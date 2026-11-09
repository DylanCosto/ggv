//! The window's side of `ggv mcp --attach` (see `ggv_io::attach` for the protocol).
//!
//! A listener thread takes loopback connections, one thread each. `state` is answered from a
//! snapshot the window refreshes every frame, so reading never waits for a frame (a minimised
//! window may not draw); commands that change the view are queued for the window, which wakes,
//! does them on its next frame and replies with its new state. What the person does (outline
//! a base, select a region, bookmark, decide a triage site) is logged as events, and `events`
//! waits for the next ones, so an agent can answer as the person points.

use std::{
    io::{BufRead, BufReader, Write},
    net::{TcpListener, TcpStream},
    collections::VecDeque,
    sync::{
        Arc, Condvar, Mutex,
        atomic::{AtomicU8, Ordering},
        mpsc::{self, Receiver, Sender},
    },
    time::{Duration, Instant},
};

use eframe::egui;
use ggv_io::attach::{Endpoint, new_token, token_matches};
use serde_json::{Value, json};

/// A request for the window, with where to send its answer.
pub struct Command {
    pub method: String,
    pub params: Value,
    pub reply: Sender<Result<Value, String>>,
    /// PENDING until either the window takes it (TAKEN) or the agent gives up (ABANDONED);
    /// exactly one wins, so a command that timed out never runs later.
    state: Arc<AtomicU8>,
}

const PENDING: u8 = 0;
const TAKEN: u8 = 1;
const ABANDONED: u8 = 2;

/// How long an agent waits for the window to take a command, and then for it to finish.
const TAKE_WAIT: Duration = Duration::from_secs(8);
const FINISH_WAIT: Duration = Duration::from_secs(60);

impl Command {
    /// Claims the command for the window; false if the agent already gave up on it, in which
    /// case it must not run.
    pub fn take(&self) -> bool {
        self.state.compare_exchange(PENDING, TAKEN, Ordering::AcqRel, Ordering::Acquire).is_ok()
    }
}

struct Snapshot {
    state: Value,
    last_used: Instant,
    focused: bool,
    /// The latest events, each with its `seq` (1, 2, ...), oldest first.
    events: VecDeque<Value>,
    seq: u64,
}

/// Events kept for agents that ask late.
const KEEP_EVENTS: usize = 200;
/// The longest an `events` call waits.
const MAX_WAIT: Duration = Duration::from_secs(300);

/// The snapshot, and a signal for agents waiting on events.
struct Shared {
    snap: Mutex<Snapshot>,
    news: Condvar,
}

pub struct Server {
    /// None for a test's server, which publishes no endpoint.
    endpoint: Option<Endpoint>,
    shared: Arc<Shared>,
    pub commands: Receiver<Command>,
}

fn new_shared() -> Arc<Shared> {
    Arc::new(Shared {
        snap: Mutex::new(Snapshot { state: json!({}), last_used: Instant::now(), focused: true, events: VecDeque::new(), seq: 0 }),
        news: Condvar::new(),
    })
}

#[cfg(test)]
impl Server {
    /// A server no agent can reach: commands go in through the sender, and `snapshot` is what
    /// an agent's `state` would read.
    pub fn detached() -> (Self, Sender<Command>) {
        let (tx, commands) = mpsc::channel();
        (Self { endpoint: None, shared: new_shared(), commands }, tx)
    }

    pub fn snapshot(&self) -> Value {
        self.shared.snap.lock().unwrap().state.clone()
    }
}

#[cfg(test)]
impl Command {
    /// A command as an agent sends it, with where its reply arrives.
    pub fn new(method: &str, params: Value) -> (Self, Receiver<Result<Value, String>>) {
        let (reply, rx) = mpsc::channel();
        (Self { method: method.to_owned(), params, reply, state: Arc::new(AtomicU8::new(PENDING)) }, rx)
    }
}

impl Server {
    /// Starts listening and publishes the endpoint file; `None` (with a message) if either fails,
    /// which leaves the app working, just not attachable.
    pub fn start(ctx: &egui::Context) -> Option<Self> {
        let listener = match TcpListener::bind(("127.0.0.1", 0)) {
            Ok(l) => l,
            Err(e) => {
                eprintln!("ggv: agents can't attach to this window: {e}");
                return None;
            }
        };
        let port = listener.local_addr().ok()?.port();
        let started = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_secs());
        let endpoint = Endpoint { pid: std::process::id(), port, token: new_token(), started };
        if let Err(e) = endpoint.publish() {
            eprintln!("ggv: agents can't attach to this window: {e}");
            return None;
        }
        let shared = new_shared();
        let (tx, commands) = mpsc::channel();
        let (token, snap, ctx) = (endpoint.token.clone(), shared.clone(), ctx.clone());
        std::thread::Builder::new()
            .name("ggv-attach".into())
            .spawn(move || {
                for stream in listener.incoming().flatten() {
                    let (token, snap, tx, ctx) = (token.clone(), snap.clone(), tx.clone(), ctx.clone());
                    let _ = std::thread::Builder::new().name("ggv-attach-conn".into()).spawn(move || serve(stream, &token, &snap, &tx, &ctx));
                }
            })
            .ok()?;
        Some(Self { endpoint: Some(endpoint), shared, commands })
    }

    /// Logs something the person did, waking agents waiting on `events`.
    pub fn push_event(&self, mut event: Value) {
        if let Ok(mut s) = self.shared.snap.lock() {
            s.seq += 1;
            event["seq"] = json!(s.seq);
            event["at_ms"] = json!(std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map_or(0, |d| d.as_millis() as u64));
            s.events.push_back(event);
            while s.events.len() > KEEP_EVENTS {
                s.events.pop_front();
            }
        }
        self.shared.news.notify_all();
    }

    /// Called every frame with the window's state.
    pub fn publish(&self, state: Value, used: bool, focused: bool) {
        if let Ok(mut s) = self.shared.snap.lock() {
            s.state = state;
            s.focused = focused;
            if used {
                s.last_used = Instant::now();
            }
        }
    }
}

impl Drop for Server {
    fn drop(&mut self) {
        if let Some(e) = &self.endpoint {
            e.unpublish();
        }
    }
}

/// The events after `since` (of `kinds`, if given), waiting up to `wait` for one; and the
/// latest event number, to pass as `since` next time.
fn events(shared: &Shared, params: &Value) -> Value {
    let kinds: Vec<String> = params.get("kinds").and_then(Value::as_array).map(|a| a.iter().filter_map(Value::as_str).map(String::from).collect()).unwrap_or_default();
    let wait = Duration::from_millis(params.get("wait_ms").and_then(Value::as_u64).unwrap_or(0)).min(MAX_WAIT);
    let deadline = Instant::now() + wait;
    let Ok(mut s) = shared.snap.lock() else { return json!({"error": "window state unavailable"}) };
    // Without `since`, only what happens from now on.
    let since = params.get("since").and_then(Value::as_u64).unwrap_or(s.seq);
    loop {
        let found: Vec<Value> = s
            .events
            .iter()
            .filter(|e| e["seq"].as_u64().is_some_and(|n| n > since))
            .filter(|e| kinds.is_empty() || e["kind"].as_str().is_some_and(|k| kinds.iter().any(|x| x == k)))
            .cloned()
            .collect();
        let now = Instant::now();
        if !found.is_empty() || now >= deadline {
            return json!({"result": {"events": found, "cursor": s.seq, "waited_ms": (wait.saturating_sub(deadline.saturating_duration_since(now))).as_millis() as u64}});
        }
        s = match shared.news.wait_timeout(s, deadline - now) {
            Ok((guard, _)) => guard,
            Err(_) => return json!({"error": "window state unavailable"}),
        };
    }
}

fn serve(stream: TcpStream, token: &str, shared: &Shared, tx: &Sender<Command>, ctx: &egui::Context) {
    let Ok(mut writer) = stream.try_clone() else { return };
    let reader = BufReader::new(stream);
    for line in reader.lines() {
        let Ok(line) = line else { return };
        if line.trim().is_empty() {
            continue;
        }
        let reply = match serde_json::from_str::<Value>(&line) {
            Err(e) => json!({"error": format!("bad request: {e}")}),
            Ok(msg) if !msg.get("token").and_then(Value::as_str).is_some_and(|t| token_matches(t, token)) => {
                // A wrong token ends the connection.
                let _ = writer.write_all(b"{\"error\":\"wrong token\"}\n");
                return;
            }
            Ok(msg) => {
                let method = msg.get("method").and_then(Value::as_str).unwrap_or("").to_owned();
                let params = msg.get("params").cloned().unwrap_or(Value::Null);
                match method.as_str() {
                    "events" => events(shared, &params),
                    "state" => match shared.snap.lock() {
                        Ok(s) => {
                            let mut v = s.state.clone();
                            v["idle_ms"] = json!(s.last_used.elapsed().as_millis() as u64);
                            v["focused"] = json!(s.focused);
                            json!({"result": v})
                        }
                        Err(_) => json!({"error": "window state unavailable"}),
                    },
                    m if crate::app::CAPABILITIES.contains(&m) => {
                        let (rtx, rrx) = mpsc::channel();
                        let state = Arc::new(AtomicU8::new(PENDING));
                        if tx.send(Command { method, params, reply: rtx, state: state.clone() }).is_err() {
                            json!({"error": "the window is closing"})
                        } else {
                            // Commands run in the app's logic step, which eframe runs even while the
                            // window is hidden or minimised; this wakes it.
                            ctx.request_repaint();
                            let answer = |r: Result<Value, String>| match r {
                                Ok(v) => json!({"result": v}),
                                Err(e) => json!({"error": e}),
                            };
                            match rrx.recv_timeout(TAKE_WAIT) {
                                Ok(r) => answer(r),
                                // Not taken yet: withdraw it, so it can't happen later unseen.
                                Err(_) if state.compare_exchange(PENDING, ABANDONED, Ordering::AcqRel, Ordering::Acquire).is_ok() => json!({
                                    "error": format!("the window didn't take the command within {} s, so it was withdrawn and not done; it's safe to retry", TAKE_WAIT.as_secs()),
                                    "completed": false,
                                }),
                                // Taken and running (a large file opening): wait for it to finish.
                                Err(_) => match rrx.recv_timeout(FINISH_WAIT) {
                                    Ok(r) => answer(r),
                                    Err(_) => json!({
                                        "error": format!("the window started the command but hasn't finished after {} s; read `state` to see whether it took effect before retrying", (TAKE_WAIT + FINISH_WAIT).as_secs()),
                                        "completed": null,
                                    }),
                                },
                            }
                        }
                    }
                    other => json!({"error": format!("unknown method `{other}`")}),
                }
            }
        };
        let mut out = reply.to_string();
        out.push('\n');
        if writer.write_all(out.as_bytes()).is_err() {
            return;
        }
    }
}
