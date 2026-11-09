//! `ggv mcp`: genome data for AI agents over the Model Context Protocol, headless.
//!
//! The protocol is newline-delimited JSON-RPC 2.0 on stdin/stdout (the MCP stdio transport);
//! only `initialize`, `tools/list`, `tools/call` and `ping` are needed, so it is implemented
//! here directly. Logging goes to stderr, never stdout.
//!
//! Input is read on its own thread, so a long call doesn't stop the server hearing the next
//! message: `ping` is answered at once, `events` waits (attached) run on threads of their own,
//! and `notifications/cancelled` stops a wait, or drops the answer of a call not yet answered.
//! Other requests are handled in order, one at a time.
//!
//! Agents don't get pictures; they get the numbers GGV draws its pictures from, computed by the
//! same engines as the desktop app and `ggv query` (so all three agree), as compact JSON with the
//! files, filters and region that produced them. See [`tools`] for what each tool returns.

pub mod desktop;
#[cfg(test)]
mod tests_tools;
pub mod tools;
mod validate;

use std::{
    collections::HashMap,
    io::{BufRead, Write},
    sync::{Arc, Mutex, PoisonError, mpsc},
    thread::JoinHandle,
};

use serde_json::{Value, json};

use desktop::Cancel;

/// Protocol versions this server speaks, newest first.
const VERSIONS: &[&str] = &["2025-06-18", "2025-03-26", "2024-11-05"];

const INSTRUCTIONS: &str = "GGV answers questions about aligned reads (BAM/CRAM), variants (VCF/BCF) and genes \
(GTF/GFF3/BED) with numbers, not pictures. Start with `info` to see the files and contigs. For a base, \
`evidence` gives allele counts by strand, read-level quality per allele, warnings for common artifacts, the VCF \
records there and what each change does to the genes (HGVS). Drill down with `reads` (reads grouped by \
allele/strand/haplotype, raw reads only if asked). For a region: `summarize`, `coverage`, `variants`, `genes`, \
`phase` (haplotypes from the reads' own alleles) and `compare` (allele-fraction differences between two files). \
Loci are like chr17:43,044,295 or chr17:43,044,295-43,125,483 (1-based, inclusive). File arguments are optional \
when the server was started with defaults (see `info`). `filter` takes an expression such as \
\"mapq >= 20 && !dup\" and applies to every count, exactly as in the GGV desktop app. When the server is attached \
to the GGV desktop app (`ggv mcp --attach`), `view` shows what the person is looking at and `navigate` moves the view; \
calls may leave out files, the base, the region and filters, which then come from the app (each answer's `desktop` \
block says what was used), so questions like \"is this variant real?\" need no coordinates.";

pub fn run(state: tools::State) -> anyhow::Result<()> {
    serve(state, std::io::BufReader::new(std::io::stdin()), std::io::stdout())
}

/// Replies go out whole, one per line, from whichever thread has one.
#[derive(Clone)]
struct Out(Arc<Mutex<Box<dyn Write + Send>>>);

impl Out {
    fn send(&self, reply: &str) -> std::io::Result<()> {
        let mut line = Vec::with_capacity(reply.len() + 1);
        line.extend_from_slice(reply.as_bytes());
        line.push(b'\n');
        let mut out = self.0.lock().unwrap_or_else(PoisonError::into_inner);
        out.write_all(&line)?;
        out.flush()
    }
}

/// Requests not yet answered, by id (as JSON text), so a cancellation can find them.
#[derive(Clone, Default)]
struct Pending(Arc<Mutex<HashMap<String, Cancel>>>);

impl Pending {
    fn add(&self, key: &str) -> Cancel {
        let c = Cancel::default();
        self.0.lock().unwrap_or_else(PoisonError::into_inner).insert(key.to_owned(), c.clone());
        c
    }

    fn cancel(&self, key: &str) {
        if let Some(c) = self.0.lock().unwrap_or_else(PoisonError::into_inner).get(key) {
            c.cancel();
        }
    }

    fn cancel_all(&self) {
        for c in self.0.lock().unwrap_or_else(PoisonError::into_inner).values() {
            c.cancel();
        }
    }

    /// Done with `key`; true if its answer should be sent (it wasn't cancelled).
    fn finish(&self, key: &str, cancel: &Cancel) -> bool {
        let mut m = self.0.lock().unwrap_or_else(PoisonError::into_inner);
        if m.get(key).is_some_and(|c| c.same(cancel)) {
            m.remove(key);
        }
        !cancel.is_cancelled()
    }
}

type Request = (String, Option<(String, Cancel)>);

/// [`run`] on any input and output.
fn serve(state: tools::State, input: impl BufRead + Send + 'static, output: impl Write + Send + 'static) -> anyhow::Result<()> {
    let out = Out(Arc::new(Mutex::new(Box::new(output))));
    let pending = Pending::default();
    // Attached: each `events` wait gets a `Desktop` of its own, so it holds nothing up.
    let waiter = state.desktop.as_ref().map(desktop::Desktop::fresh);
    let (tx, rx) = mpsc::channel::<Request>();
    let reader = {
        let (out, pending) = (out.clone(), pending.clone());
        std::thread::Builder::new().name("ggv-mcp-input".into()).spawn(move || read(input, &out, &pending, waiter.as_ref(), &tx))?
    };
    let mut server = Server { state };
    for (line, req) in rx {
        if let Some((key, cancel)) = &req
            && cancel.is_cancelled()
        {
            pending.finish(key, cancel);
            continue;
        }
        let reply = server.handle(&line);
        let send = req.as_ref().is_none_or(|(key, cancel)| pending.finish(key, cancel));
        if let Some(reply) = reply
            && send
            && let Err(e) = out.send(&reply)
        {
            // The client has gone: stop the waits too. (The input thread is left blocked on a
            // read; returning ends the process.)
            pending.cancel_all();
            return Err(e.into());
        }
    }
    reader.join().map_err(|_| anyhow::anyhow!("the input thread failed"))?
}

/// Reads messages until the input ends, answering `ping` and cancellations itself, starting
/// `events` waits, and passing the rest on to be handled in order. Waits still going when the
/// input ends are stopped, unanswered.
fn read(input: impl BufRead, out: &Out, pending: &Pending, waiter: Option<&desktop::Desktop>, tx: &mpsc::Sender<Request>) -> anyhow::Result<()> {
    let mut waits: Vec<(Cancel, JoinHandle<()>)> = Vec::new();
    let mut result = Ok(());
    let mut input = input;
    let mut bytes = Vec::new();
    loop {
        bytes.clear();
        match input.read_until(b'\n', &mut bytes) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => {
                result = Err(e.into());
                break;
            }
        }
        // One line that isn't UTF-8 is a parse error, not the end of the session.
        let line = match std::str::from_utf8(&bytes) {
            Ok(l) => l.trim_end_matches(['\n', '\r']).to_string(),
            Err(e) => {
                if out.send(&error(&Value::Null, -32700, format!("parse error: not UTF-8: {e}")).to_string()).is_err() {
                    break;
                }
                continue;
            }
        };
        if line.trim().is_empty() {
            continue;
        }
        waits.retain(|(_, h)| !h.is_finished());
        let msg: Value = serde_json::from_str(&line).unwrap_or(Value::Null);
        let method = msg.get("method").and_then(Value::as_str).unwrap_or("");
        let id = msg.get("id");
        match (method, id) {
            ("ping", Some(id)) => {
                if out.send(&json!({"jsonrpc": "2.0", "id": id, "result": {}}).to_string()).is_err() {
                    break;
                }
                continue;
            }
            ("notifications/cancelled", None) => {
                if let Some(r) = msg.pointer("/params/requestId") {
                    pending.cancel(&r.to_string());
                }
                continue;
            }
            _ => {}
        }
        let wait_on = waiter.filter(|_| method == "tools/call" && msg.pointer("/params/name") == Some(&json!("events")));
        if let (Some(waiter), Some(id)) = (wait_on, id) {
            let key = id.to_string();
            let cancel = pending.add(&key);
            let mut desktop = waiter.fresh();
            let (id, out, pending, c) = (id.clone(), out.clone(), pending.clone(), cancel.clone());
            let args = msg.pointer("/params/arguments").cloned().unwrap_or_else(|| json!({}));
            let h = std::thread::Builder::new().name("ggv-mcp-events".into()).spawn(move || {
                let t = std::time::Instant::now();
                let window = args.get("window").and_then(Value::as_u64).map(|w| w as u32);
                let r = desktop.events_tool(window, &args, &c).map(tools::Output::Json).map_err(|e| anyhow::anyhow!(e));
                eprintln!("ggv mcp: events in {:.1} ms", t.elapsed().as_secs_f64() * 1000.0);
                // A wait stopped because the input ended isn't answered; one that finished is.
                if pending.finish(&key, &c) && !c.interrupted() {
                    let _ = out.send(&json!({"jsonrpc": "2.0", "id": id, "result": tool_result(r)}).to_string());
                }
            })?;
            waits.push((cancel, h));
            continue;
        }
        let req = id.map(|id| {
            let key = id.to_string();
            let c = pending.add(&key);
            (key, c)
        });
        if tx.send((line, req)).is_err() {
            break;
        }
    }
    for (c, h) in waits {
        c.close();
        let _ = h.join();
    }
    result
}

/// A tool's outcome as a `tools/call` result.
fn tool_result(r: anyhow::Result<tools::Output>) -> Value {
    match r {
        Ok(tools::Output::Json(v)) => json!({"content": [{"type": "text", "text": v.to_string()}], "isError": false}),
        Ok(tools::Output::Image { png, note }) => json!({
            "content": [
                {"type": "image", "data": ggv_render::report::base64(&png), "mimeType": "image/png"},
                {"type": "text", "text": note},
            ],
            "isError": false
        }),
        // Tool failures go back to the model as results it can read and act on.
        Err(e) => json!({"content": [{"type": "text", "text": format!("error: {e:#}")}], "isError": true}),
    }
}

pub struct Server {
    pub state: tools::State,
}

fn error(id: &Value, code: i64, message: impl Into<String>) -> Value {
    json!({"jsonrpc": "2.0", "id": id, "error": {"code": code, "message": message.into()}})
}

impl Server {
    /// Handles one message; returns the reply (none for notifications).
    pub fn handle(&mut self, line: &str) -> Option<String> {
        let msg: Value = match serde_json::from_str(line) {
            Ok(v) => v,
            Err(e) => return Some(error(&Value::Null, -32700, format!("parse error: {e}")).to_string()),
        };
        // A batch (older protocol versions): answer each.
        if let Value::Array(items) = msg {
            if items.is_empty() {
                return Some(error(&Value::Null, -32600, "invalid request: an empty batch").to_string());
            }
            let replies: Vec<Value> = items.iter().filter_map(|m| self.one(m)).collect();
            return (!replies.is_empty()).then(|| Value::Array(replies).to_string());
        }
        self.one(&msg).map(|v| v.to_string())
    }

    fn one(&mut self, msg: &Value) -> Option<Value> {
        let Some(obj) = msg.as_object() else {
            return Some(error(&Value::Null, -32600, "invalid request: not a JSON-RPC object"));
        };
        let id = obj.get("id").cloned();
        let Some(method) = obj.get("method").and_then(Value::as_str) else {
            // A response to us (we send no requests) needs no reply; anything else is invalid.
            if obj.contains_key("result") || obj.contains_key("error") {
                return None;
            }
            return Some(error(&id.unwrap_or(Value::Null), -32600, "invalid request: no method"));
        };
        let params = obj.get("params").cloned().unwrap_or(Value::Null);
        // Notifications (no id) get no reply.
        let id = id?;
        let result = match method {
            "initialize" => {
                let asked = params.get("protocolVersion").and_then(Value::as_str).unwrap_or(VERSIONS[0]);
                let version = VERSIONS.iter().find(|v| **v == asked).copied().unwrap_or(VERSIONS[0]);
                json!({
                    "protocolVersion": version,
                    "capabilities": {"tools": {"listChanged": false}},
                    "serverInfo": {"name": "ggv", "title": "GGV genome viewer", "version": concat!(env!("CARGO_PKG_VERSION"), "+", env!("GGV_COMMIT"))},
                    "instructions": INSTRUCTIONS,
                })
            }
            "ping" => json!({}),
            "tools/list" => json!({"tools": self.state.list()}),
            "tools/call" => {
                let Some(name) = params.get("name").and_then(Value::as_str) else {
                    return Some(error(&id, -32602, "tools/call needs a tool name"));
                };
                let args = params.get("arguments").cloned().unwrap_or_else(|| json!({}));
                let tools = self.state.list();
                let Some(tool) = tools.iter().find(|t| t["name"] == name) else {
                    return Some(error(&id, -32602, format!("unknown tool: {name}")));
                };
                // Supplied arguments must fit the published schema: a wrong type or an
                // out-of-range value is refused, never replaced by a default.
                if let Err(e) = validate::check(&tool["inputSchema"], &args) {
                    return Some(json!({"jsonrpc": "2.0", "id": id, "result": {"content": [{"type": "text", "text": format!("error: invalid arguments: {e}")}], "isError": true}}));
                }
                let t = std::time::Instant::now();
                let result = tool_result(self.state.call(name, &args));
                eprintln!("ggv mcp: {name} in {:.1} ms", t.elapsed().as_secs_f64() * 1000.0);
                result
            }
            "resources/list" => json!({"resources": []}),
            "prompts/list" => json!({"prompts": []}),
            _ => return Some(error(&id, -32601, format!("method not found: {method}"))),
        };
        Some(json!({"jsonrpc": "2.0", "id": id, "result": result}))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn server() -> Server {
        Server { state: tools::State::new(tools::Defaults::default()) }
    }

    fn call(s: &mut Server, msg: Value) -> Value {
        serde_json::from_str(&s.handle(&msg.to_string()).expect("a reply")).unwrap()
    }

    #[test]
    fn handshake_and_listing() {
        let mut s = server();
        let r = call(&mut s, json!({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": "2025-03-26", "capabilities": {}, "clientInfo": {"name": "t", "version": "0"}}}));
        assert_eq!(r["result"]["protocolVersion"], "2025-03-26");
        assert_eq!(r["result"]["serverInfo"]["name"], "ggv");
        assert!(s.handle(r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#).is_none(), "notifications get no reply");
        let r = call(&mut s, json!({"jsonrpc": "2.0", "id": "x", "method": "initialize", "params": {"protocolVersion": "1999-01-01"}}));
        assert_eq!(r["result"]["protocolVersion"], VERSIONS[0], "unknown versions get ours");
        let tools = call(&mut s, json!({"jsonrpc": "2.0", "id": 2, "method": "tools/list"}));
        let names: Vec<&str> = tools["result"]["tools"].as_array().unwrap().iter().map(|t| t["name"].as_str().unwrap()).collect();
        for want in ["info", "evidence", "reads", "summarize", "coverage", "variants", "genes", "consequence", "phase", "compare", "render", "open"] {
            assert!(names.contains(&want), "{want} missing from {names:?}");
        }
        for t in tools["result"]["tools"].as_array().unwrap() {
            assert_eq!(t["inputSchema"]["type"], "object", "{}", t["name"]);
            assert!(t["description"].as_str().unwrap().len() > 20);
        }
    }

    #[test]
    fn errors_are_reported_not_fatal() {
        let mut s = server();
        let r = call(&mut s, json!({"jsonrpc": "2.0", "id": 3, "method": "nope"}));
        assert_eq!(r["error"]["code"], -32601);
        let r = call(&mut s, json!({"jsonrpc": "2.0", "id": 4, "method": "tools/call", "params": {"name": "nope"}}));
        assert_eq!(r["error"]["code"], -32602);
        // A tool that fails (no files) answers with isError and a message the model can use.
        let r = call(&mut s, json!({"jsonrpc": "2.0", "id": 5, "method": "tools/call", "params": {"name": "evidence", "arguments": {"locus": "chr1:100"}}}));
        assert_eq!(r["result"]["isError"], true);
        assert!(r["result"]["content"][0]["text"].as_str().unwrap().contains("no alignment files"));
        assert_eq!(serde_json::from_str::<Value>(&s.handle("{not json").unwrap()).unwrap()["error"]["code"], -32700);
        assert_eq!(call(&mut s, json!({"jsonrpc": "2.0", "id": 6, "method": "ping"}))["result"], json!({}));
        // Not a request: an empty batch, a bare value, no method.
        for (msg, id) in [("[]", Value::Null), ("5", Value::Null), (r#"{"jsonrpc":"2.0","id":7}"#, json!(7))] {
            let r: Value = serde_json::from_str(&s.handle(msg).expect("a reply")).unwrap();
            assert_eq!((&r["error"]["code"], &r["id"]), (&json!(-32600), &id), "{msg}");
        }
        let r: Value = serde_json::from_str(&s.handle("[5]").unwrap()).unwrap();
        assert_eq!(r[0]["error"]["code"], -32600);
        assert!(s.handle(r#"{"jsonrpc":"2.0","id":8,"result":{}}"#).is_none(), "a response needs no reply");
    }

    #[test]
    fn consequence_finds_transcripts_the_deletion_reaches() {
        // chrQA:100-102 is AAC; the transcript starts at 101. AAC>A at 100 deletes 101-102.
        let dir = std::env::temp_dir().join(format!("ggv-mcp-consequence-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let seq = format!("{}AAC{}", "T".repeat(99), "G".repeat(98));
        std::fs::write(dir.join("ref.fa"), format!(">chrQA\n{seq}\n")).unwrap();
        std::fs::write(dir.join("ref.fa.fai"), "chrQA\t200\t7\t200\t201\n").unwrap();
        let gff = "##gff-version 3\nchrQA\tqa\tgene\t101\t180\t.\t+\t.\tID=g;Name=QA\nchrQA\tqa\tmRNA\t101\t180\t.\t+\t.\tID=t;Parent=g\nchrQA\tqa\texon\t101\t180\t.\t+\t.\tParent=t\n";
        std::fs::write(dir.join("genes.gff3"), gff).unwrap();
        let mut s = server();
        let args = json!({"locus": "chrQA:100", "ref": "AAC", "alt": "A", "reference": dir.join("ref.fa"), "genes": [dir.join("genes.gff3")]});
        let r = call(&mut s, json!({"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": "consequence", "arguments": args}}));
        let v: Value = serde_json::from_str(r["result"]["content"][0]["text"].as_str().unwrap()).unwrap();
        assert_eq!(v["consequences"][0]["gene"], "QA", "{v}");
        assert_eq!(v["consequences"][0]["hgvs_c"], "n.1_2del");
        std::fs::remove_dir_all(&dir).unwrap();
    }
}

#[cfg(test)]
mod transport_tests {
    use std::{
        io::{BufReader, Read},
        net::TcpListener,
        sync::mpsc::{Receiver, Sender, channel},
        time::{Duration, Instant},
    };

    use super::*;

    /// Input fed line by line from the test; dropping the sender ends it.
    struct Input(Receiver<Vec<u8>>, Vec<u8>);

    impl Read for Input {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            if self.1.is_empty() {
                match self.0.recv() {
                    Ok(v) => self.1 = v,
                    Err(_) => return Ok(0),
                }
            }
            let n = buf.len().min(self.1.len());
            buf[..n].copy_from_slice(&self.1[..n]);
            self.1.drain(..n);
            Ok(n)
        }
    }

    /// Each reply arrives whole.
    struct Output(Sender<Vec<u8>>);

    impl Write for Output {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            let _ = self.0.send(buf.to_vec());
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    /// A stand-in GGV window: answers `state`, and holds `events` for its `wait_ms` unless the
    /// connection is closed first, which it reports.
    fn fake_window(dir: &std::path::Path, pid: u32) -> Receiver<()> {
        let listener = TcpListener::bind(("127.0.0.1", 0)).unwrap();
        let ep = ggv_io::attach::Endpoint { pid, port: listener.local_addr().unwrap().port(), token: "t".into(), started: 1 };
        ep.publish_in(dir).unwrap();
        let (tx, rx) = channel();
        std::thread::spawn(move || {
            for stream in listener.incoming().flatten() {
                let tx = tx.clone();
                std::thread::spawn(move || {
                    let mut w = stream.try_clone().unwrap();
                    let mut r = BufReader::new(stream);
                    let mut line = String::new();
                    while r.read_line(&mut line).unwrap_or(0) > 0 {
                        let req: Value = serde_json::from_str(&line).unwrap();
                        line.clear();
                        if req["method"] == "events" {
                            let wait = req["params"]["wait_ms"].as_u64().unwrap();
                            if wait > 0 {
                                r.get_ref().set_read_timeout(Some(Duration::from_millis(wait))).unwrap();
                                match r.read_line(&mut line) {
                                    Err(e) if matches!(e.kind(), std::io::ErrorKind::WouldBlock | std::io::ErrorKind::TimedOut) => {}
                                    _ => {
                                        let _ = tx.send(());
                                        return;
                                    }
                                }
                            }
                            let _ = writeln!(w, "{}", json!({"result": {"events": [], "cursor": 0}}));
                        } else {
                            let _ = writeln!(w, "{}", json!({"result": {"panes": [], "focused": true}}));
                        }
                    }
                });
            }
        });
        rx
    }

    #[test]
    fn a_line_that_is_not_utf8_is_a_parse_error_not_the_end() {
        let (to_server, input) = channel::<Vec<u8>>();
        let (output, replies) = channel::<Vec<u8>>();
        let server = std::thread::spawn(move || serve(tools::State::new(tools::Defaults::default()), BufReader::new(Input(input, Vec::new())), Output(output)));
        let reply = || -> Value { serde_json::from_slice(&replies.recv_timeout(Duration::from_secs(10)).expect("a reply in time")).unwrap() };
        to_server.send(b"{\"jsonrpc\":\"2.0\",\"id\":1,\"method\":\"tools/call\",\"params\":{\"name\":\"info\",\"arguments\":{\"files\":[\"C:/x\xe9.bam\"]}}}\n".to_vec()).unwrap();
        let r = reply();
        assert_eq!((&r["error"]["code"], &r["id"]), (&json!(-32700), &Value::Null), "{r}");
        to_server.send(b"{\"jsonrpc\":\"2.0\",\"id\":2,\"method\":\"ping\"}\r\n".to_vec()).unwrap();
        assert_eq!(reply()["id"], 2);
        drop(to_server);
        server.join().unwrap().unwrap();
    }

    #[test]
    fn ping_answers_during_a_wait_and_cancelled_waits_stay_silent() {
        let dir = std::env::temp_dir().join(format!("ggv-mcp-transport-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let pid = 4_000_000_002;
        let interrupted = fake_window(&dir, pid);
        let mut state = tools::State::new(tools::Defaults::default());
        state.desktop = Some(desktop::Desktop::in_dir(Some(pid), dir.clone()));
        let (to_server, input) = channel::<Vec<u8>>();
        let (output, replies) = channel::<Vec<u8>>();
        let server = std::thread::spawn(move || serve(state, BufReader::new(Input(input, Vec::new())), Output(output)));
        let send = |v: Value| to_server.send(format!("{v}\n").into_bytes()).unwrap();
        let reply = |within: Duration| -> Value { serde_json::from_slice(&replies.recv_timeout(within).expect("a reply in time")).unwrap() };

        send(json!({"jsonrpc": "2.0", "id": 0, "method": "initialize", "params": {"protocolVersion": "2025-06-18"}}));
        assert_eq!(reply(Duration::from_secs(10))["id"], 0);

        // A long wait doesn't hold up a ping, or another call.
        send(json!({"jsonrpc": "2.0", "id": 1, "method": "tools/call", "params": {"name": "events", "arguments": {"wait_s": 30}}}));
        std::thread::sleep(Duration::from_millis(200));
        let t = Instant::now();
        send(json!({"jsonrpc": "2.0", "id": 2, "method": "ping"}));
        assert_eq!(reply(Duration::from_secs(5))["id"], 2);
        assert!(t.elapsed() < Duration::from_secs(2), "ping took {:?}", t.elapsed());
        send(json!({"jsonrpc": "2.0", "id": 3, "method": "tools/list"}));
        assert_eq!(reply(Duration::from_secs(5))["id"], 3);

        // Cancelling ends the wait at once (the window sees the connection close) and sends no
        // result for it.
        send(json!({"jsonrpc": "2.0", "method": "notifications/cancelled", "params": {"requestId": 1, "reason": "test"}}));
        interrupted.recv_timeout(Duration::from_secs(5)).expect("the wait was interrupted");
        send(json!({"jsonrpc": "2.0", "id": 4, "method": "ping"}));
        assert_eq!(reply(Duration::from_secs(5))["id"], 4, "no answer for the cancelled call");

        // A wait that isn't cancelled still answers.
        send(json!({"jsonrpc": "2.0", "id": 5, "method": "tools/call", "params": {"name": "events", "arguments": {"wait_s": 0.2}}}));
        let r = reply(Duration::from_secs(10));
        assert_eq!(r["id"], 5);
        assert_eq!(r["result"]["isError"], false, "{r}");

        // Closing the input ends a wait still going, and the server, without answering it.
        send(json!({"jsonrpc": "2.0", "id": 6, "method": "tools/call", "params": {"name": "events", "arguments": {"wait_s": 30}}}));
        std::thread::sleep(Duration::from_millis(200));
        drop(to_server);
        let t = Instant::now();
        server.join().unwrap().unwrap();
        assert!(t.elapsed() < Duration::from_secs(5));
        interrupted.recv_timeout(Duration::from_secs(5)).expect("the wait was interrupted");
        assert!(replies.try_recv().is_err(), "nothing more was sent");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
