//! Attaching to a running GGV desktop window, for `ggv mcp --attach`.
//!
//! Each desktop window listens on a loopback TCP port (127.0.0.1 only) and writes an endpoint
//! file, `<config>/ggv/attach/<pid>.json`, with its port and a random token; the file lives in the
//! user's own config folder, and every request must carry the token, so other users and other
//! machines can't talk to the window. The window removes the file when it closes; a file whose
//! window has gone (a crash) is noticed when connecting fails, and removed.
//!
//! The protocol is newline-delimited JSON: a request `{"token", "method", "params"}`, a reply
//! `{"result": ...}` or `{"error": "..."}`, one per line, on a connection kept open.
//!
//! Methods: `state` (files, panes, view, outlined base, selection, filters; answered without
//! waiting for the window to draw), `goto` / `select` / `open` (done by the window on its next
//! frame; the reply carries the new state).

use std::{
    io::{BufRead, BufReader, ErrorKind, Write},
    net::{SocketAddr, TcpStream},
    path::{Path, PathBuf},
    time::Duration,
};

use serde::{Deserialize, Serialize};

/// Where a window can be reached.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Endpoint {
    pub pid: u32,
    pub port: u16,
    pub token: String,
    /// When the window started (Unix seconds).
    pub started: u64,
}

/// The per-user config folder GGV uses (`%APPDATA%\ggv`, `$XDG_CONFIG_HOME/ggv`, `~/.config/ggv`).
pub fn config_dir() -> Option<PathBuf> {
    crate::folders::config_dir()
}

/// The folder of endpoint files (`GGV_ATTACH_DIR` overrides it, for tests).
pub fn endpoint_dir() -> Option<PathBuf> {
    match std::env::var_os("GGV_ATTACH_DIR") {
        Some(d) => Some(PathBuf::from(d)),
        None => config_dir().map(|d| d.join("attach")),
    }
}

/// A random token: 128 bits from the standard library's OS-seeded hasher keys.
pub fn new_token() -> String {
    use std::hash::{BuildHasher, Hasher};
    let nanos = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_nanos());
    (0..2u64)
        .map(|i| {
            let mut h = std::collections::hash_map::RandomState::new().build_hasher();
            h.write_u64(i);
            h.write_u128(nanos);
            h.write_u32(std::process::id());
            format!("{:016x}", h.finish())
        })
        .collect()
}

impl Endpoint {
    pub fn path(&self) -> Option<PathBuf> {
        endpoint_dir().map(|d| d.join(format!("{}.json", self.pid)))
    }

    /// Writes the endpoint file, readable by this user only: on Unix the folder is made 0700
    /// and the file 0600 whatever the umask; on Windows both sit in the user's own profile
    /// (`%APPDATA%`), which only that user (and administrators) can read.
    pub fn publish(&self) -> std::io::Result<PathBuf> {
        let dir = endpoint_dir().ok_or_else(|| std::io::Error::other("no config folder"))?;
        self.publish_in(&dir)
    }

    /// [`publish`](Self::publish) into `dir`.
    pub fn publish_in(&self, dir: &Path) -> std::io::Result<PathBuf> {
        std::fs::create_dir_all(dir)?;
        let path = dir.join(format!("{}.json", self.pid));
        // Written whole under a temporary name, then renamed, so readers never see half a file.
        // `create_new` won't follow a link left at the temporary name.
        let tmp = dir.join(format!("{}.json.tmp", self.pid));
        let _ = std::fs::remove_file(&tmp);
        let mut open = std::fs::OpenOptions::new();
        open.write(true).create_new(true);
        #[cfg(unix)]
        std::os::unix::fs::OpenOptionsExt::mode(&mut open, 0o600);
        let mut f = open.open(&tmp)?;
        let written = (|| {
            #[cfg(unix)]
            restrict(dir, &f)?;
            f.write_all(&serde_json::to_vec(self).map_err(std::io::Error::other)?)?;
            f.sync_all()?;
            drop(f);
            std::fs::rename(&tmp, &path)
        })();
        if written.is_err() {
            let _ = std::fs::remove_file(&tmp);
        }
        written.map(|()| path)
    }

    pub fn unpublish(&self) {
        if let Some(p) = self.path() {
            let _ = std::fs::remove_file(p);
        }
    }
}

/// Makes the endpoint folder 0700 and the new file `f` 0600, whatever the umask; a folder owned
/// by someone else is refused rather than trusted.
#[cfg(unix)]
fn restrict(dir: &Path, f: &std::fs::File) -> std::io::Result<()> {
    use std::os::unix::fs::{MetadataExt, PermissionsExt};
    f.set_permissions(std::fs::Permissions::from_mode(0o600))?;
    // The new file belongs to this process's user.
    let me = f.metadata()?.uid();
    let d = std::fs::metadata(dir)?;
    if d.uid() != me {
        return Err(std::io::Error::new(ErrorKind::PermissionDenied, format!("{} belongs to another user", dir.display())));
    }
    if d.mode() & 0o777 != 0o700 {
        std::fs::set_permissions(dir, std::fs::Permissions::from_mode(0o700))?;
    }
    Ok(())
}

/// Every endpoint file, newest window first. Unreadable files are skipped.
pub fn endpoints() -> Vec<Endpoint> {
    endpoint_dir().map_or_else(Vec::new, |d| endpoints_in(&d))
}

/// [`endpoints`] in `dir`.
pub fn endpoints_in(dir: &Path) -> Vec<Endpoint> {
    let Ok(entries) = std::fs::read_dir(dir) else { return Vec::new() };
    let mut v: Vec<Endpoint> = entries
        .filter_map(Result::ok)
        .filter(|e| e.path().extension().is_some_and(|x| x == "json"))
        .filter_map(|e| serde_json::from_slice(&std::fs::read(e.path()).ok()?).ok())
        .collect();
    v.sort_by(|a, b| b.started.cmp(&a.started).then(b.pid.cmp(&a.pid)));
    v
}

/// The longest a window takes to answer a command (8 s to take it, 60 s to finish), plus slack.
pub const REPLY_WAIT: Duration = Duration::from_secs(75);

/// A connection to one window.
pub struct Client {
    pub endpoint: Endpoint,
    reader: BufReader<TcpStream>,
    writer: TcpStream,
}

impl Client {
    /// Connects; a window that has gone has its endpoint file removed.
    pub fn connect(endpoint: &Endpoint) -> std::io::Result<Self> {
        let addr = SocketAddr::from(([127, 0, 0, 1], endpoint.port));
        let stream = match TcpStream::connect_timeout(&addr, Duration::from_millis(500)) {
            Ok(s) => s,
            Err(e) => {
                endpoint.unpublish();
                return Err(e);
            }
        };
        stream.set_nodelay(true)?;
        // The window answers within REPLY_WAIT: it withdraws a command it couldn't take in 8 s,
        // and gives one it started (a big file opening) up to a minute more.
        stream.set_read_timeout(Some(REPLY_WAIT))?;
        Ok(Self { endpoint: endpoint.clone(), reader: BufReader::new(stream.try_clone()?), writer: stream })
    }

    /// A handle to the connection that another thread can shut down, which ends a request
    /// waiting on it at once (with an error).
    pub fn interrupter(&self) -> std::io::Result<TcpStream> {
        self.writer.try_clone()
    }

    /// [`request`](Self::request) that may wait up to `timeout` for the reply (the usual limit
    /// is [`REPLY_WAIT`]), for requests that wait on the person.
    pub fn request_waiting(&mut self, method: &str, params: serde_json::Value, timeout: Duration) -> Result<serde_json::Value, RequestError> {
        let set = |c: &Self, t| c.writer.set_read_timeout(Some(t)).map_err(|e| RequestError::stale(format!("window {}: {e}", c.endpoint.pid)));
        set(self, timeout)?;
        let r = self.request(method, params);
        let _ = set(self, REPLY_WAIT);
        r
    }

    pub fn request(&mut self, method: &str, params: serde_json::Value) -> Result<serde_json::Value, RequestError> {
        let pid = self.endpoint.pid;
        let msg = serde_json::json!({"token": self.endpoint.token, "method": method, "params": params});
        let mut line = msg.to_string();
        line.push('\n');
        self.writer.write_all(line.as_bytes()).map_err(|e| RequestError::stale(format!("window {pid}: {e}")))?;
        let mut reply = String::new();
        match self.reader.read_line(&mut reply) {
            Ok(0) => return Err(RequestError::stale(format!("window {pid} closed the connection"))),
            Ok(_) => {}
            // A reset connection was dead before the window read the request; a timeout wasn't,
            // and the window may still carry it out (its late reply would then answer the next
            // request on this connection).
            Err(e) if matches!(e.kind(), ErrorKind::ConnectionReset | ErrorKind::ConnectionAborted | ErrorKind::BrokenPipe) => {
                return Err(RequestError::stale(format!("window {pid}: {e}")));
            }
            Err(e) => return Err(RequestError::lost(format!("window {pid}: {e}"))),
        }
        let v: serde_json::Value = serde_json::from_str(&reply).map_err(|e| RequestError::lost(format!("window {pid}: bad reply: {e}")))?;
        match v.get("error").and_then(serde_json::Value::as_str) {
            Some(e) => Err(RequestError::window(e.to_owned())),
            None => Ok(v.get("result").cloned().unwrap_or(serde_json::Value::Null)),
        }
    }
}

/// Why a request failed, and whether it's safe to send again.
#[derive(Debug)]
pub struct RequestError {
    pub message: String,
    /// The connection was dead, so the window never saw the request: sending it again on a
    /// fresh connection can't do it twice. Otherwise the window has it (and may still act on it).
    pub stale: bool,
    /// The connection can't be used again either: a reply timed out (and may still come) or
    /// came garbled, so replies on it no longer line up with requests.
    pub lost: bool,
}

impl RequestError {
    fn stale(message: String) -> Self {
        Self { message, stale: true, lost: true }
    }

    fn lost(message: String) -> Self {
        Self { message, stale: false, lost: true }
    }

    fn window(message: String) -> Self {
        Self { message, stale: false, lost: false }
    }
}

impl std::fmt::Display for RequestError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

/// Compares tokens without stopping at the first difference.
pub fn token_matches(given: &str, expected: &str) -> bool {
    given.len() == expected.len() && given.bytes().zip(expected.bytes()).fold(0u8, |acc, (a, b)| acc | (a ^ b)) == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokens_are_random_and_compared_whole() {
        let (a, b) = (new_token(), new_token());
        assert_eq!(a.len(), 32);
        assert_ne!(a, b);
        assert!(token_matches(&a, &a));
        assert!(!token_matches(&a, &b));
        assert!(!token_matches(&a, &a[..31]));
    }

    #[test]
    fn endpoint_files_are_private_and_whole() {
        let dir = std::env::temp_dir().join(format!("ggv-attach-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        // A folder left world-readable (an older GGV, or umask 022) is tightened.
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&dir, std::fs::Permissions::from_mode(0o755)).unwrap();
        }
        let ep = Endpoint { pid: 4_000_000_001, port: 1234, token: new_token(), started: 1 };
        let path = ep.publish_in(&dir).unwrap();
        assert_eq!(endpoints_in(&dir), vec![ep.clone()]);
        // Published again (same pid): replaced, no temporary file left behind.
        ep.publish_in(&dir).unwrap();
        let names: Vec<_> = std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().file_name()).collect();
        assert_eq!(names, vec![std::ffi::OsString::from("4000000001.json")]);
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            assert_eq!(std::fs::metadata(&dir).unwrap().permissions().mode() & 0o777, 0o700);
            assert_eq!(std::fs::metadata(&path).unwrap().permissions().mode() & 0o777, 0o600);
        }
        let _ = path;
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
