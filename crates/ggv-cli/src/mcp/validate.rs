//! Checks a tool call's arguments against the tool's published input schema before the call
//! runs, so a supplied value of the wrong type or out of range is refused (naming the argument)
//! instead of being dropped for a default. Covers the parts of JSON Schema the tools use: type,
//! enum, minimum/maximum, items, maxItems, properties, additionalProperties, anyOf and required.
//! `null` counts as left out. Kept for compatibility: a single string where a list of strings is
//! expected is taken as a one-item list.

use serde_json::Value;

/// `Err` names the first argument that doesn't fit.
pub fn check(schema: &Value, args: &Value) -> Result<(), String> {
    let Some(obj) = args.as_object() else {
        return if args.is_null() { Ok(()) } else { Err("arguments must be an object".into()) };
    };
    for r in schema["required"].as_array().into_iter().flatten().filter_map(Value::as_str) {
        if obj.get(r).is_none_or(Value::is_null) {
            return Err(format!("`{r}` is required"));
        }
    }
    let props = &schema["properties"];
    for (k, v) in obj {
        if let Some(s) = props.get(k) {
            value(s, v, k)?;
        }
    }
    Ok(())
}

fn type_name(v: &Value) -> &'static str {
    match v {
        Value::Null => "null",
        Value::Bool(_) => "a boolean",
        Value::Number(n) if n.is_u64() || n.is_i64() => "an integer",
        Value::Number(_) => "a number",
        Value::String(_) => "a string",
        Value::Array(_) => "a list",
        Value::Object(_) => "an object",
    }
}

fn is_type(t: &str, v: &Value) -> bool {
    match t {
        "string" => v.is_string(),
        "integer" => v.is_u64() || v.is_i64(),
        "number" => v.is_number(),
        "boolean" => v.is_boolean(),
        "array" => v.is_array(),
        "object" => v.is_object(),
        "null" => v.is_null(),
        _ => true,
    }
}

fn value(s: &Value, v: &Value, at: &str) -> Result<(), String> {
    if v.is_null() {
        return Ok(());
    }
    if let Some(any) = s["anyOf"].as_array() {
        // The first branch's error is the most telling one.
        let mut first = None;
        for b in any {
            match value(b, v, at) {
                Ok(()) => {
                    first = None;
                    break;
                }
                Err(e) => {
                    first.get_or_insert(e);
                }
            }
        }
        if let Some(e) = first {
            return Err(e);
        }
    }
    let types: Vec<&str> = match &s["type"] {
        Value::String(t) => vec![t.as_str()],
        Value::Array(a) => a.iter().filter_map(Value::as_str).collect(),
        _ => Vec::new(),
    };
    let string_list = types == ["array"] && s["items"]["type"] == "string" && v.is_string();
    if !types.is_empty() && !string_list && !types.iter().any(|t| is_type(t, v)) {
        let want = types.iter().map(|t| if *t == "array" { "list" } else { t }).collect::<Vec<_>>().join(" or ");
        return Err(format!("`{at}` must be {} {want}, not {}", if want.starts_with(['a', 'i', 'o']) { "an" } else { "a" }, type_name(v)));
    }
    if let Some(e) = s["enum"].as_array()
        && !e.contains(v)
    {
        let opts: Vec<String> = e.iter().map(|x| x.as_str().map_or_else(|| x.to_string(), String::from)).collect();
        return Err(format!("`{at}` must be one of {}, not {v}", opts.join(", ")));
    }
    if let Some(x) = v.as_f64() {
        if let Some(min) = s["minimum"].as_f64()
            && x < min
        {
            return Err(format!("`{at}` must be at least {}, not {v}", s["minimum"]));
        }
        if let Some(max) = s["maximum"].as_f64()
            && x > max
        {
            return Err(format!("`{at}` must be at most {}, not {v}", s["maximum"]));
        }
    }
    if let Some(a) = v.as_array() {
        if let Some(max) = s["maxItems"].as_u64()
            && a.len() as u64 > max
        {
            return Err(format!("`{at}` takes at most {max} items, not {}", a.len()));
        }
        if s["items"].is_object() {
            for (i, x) in a.iter().enumerate() {
                if x.is_null() {
                    return Err(format!("`{at}[{i}]` is null"));
                }
                value(&s["items"], x, &format!("{at}[{i}]"))?;
            }
        }
    }
    if let Some(o) = v.as_object() {
        for (k, x) in o {
            let sub = s["properties"].get(k).or_else(|| s.get("additionalProperties").filter(|a| a.is_object()));
            if let Some(sub) = sub {
                value(sub, x, &format!("{at}.{k}"))?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn types_ranges_and_lists() {
        let s = json!({"type": "object", "required": ["region"], "properties": {
            "region": {"type": "string"},
            "min_mapq": {"type": "integer", "minimum": 0, "maximum": 255},
            "threshold": {"type": "number", "minimum": 0.01, "maximum": 1.0},
            "bams": {"type": "array", "items": {"type": "string"}},
            "mode": {"type": "string", "enum": ["a", "b"]},
            "track": {"type": ["integer", "string"]},
        }});
        assert!(check(&s, &json!({"region": "chr1"})).is_ok());
        assert!(check(&s, &json!({})).unwrap_err().contains("`region` is required"));
        assert!(check(&s, &json!({"region": "x", "min_mapq": "255"})).unwrap_err().contains("`min_mapq` must be an integer"));
        assert!(check(&s, &json!({"region": "x", "min_mapq": 256})).unwrap_err().contains("at most 255"));
        assert!(check(&s, &json!({"region": "x", "min_mapq": 2.5})).is_err());
        assert!(check(&s, &json!({"region": "x", "threshold": -0.01})).unwrap_err().contains("`threshold` must be at least"));
        assert!(check(&s, &json!({"region": "x", "threshold": 1})).is_ok());
        assert!(check(&s, &json!({"region": "x", "bams": [true]})).unwrap_err().contains("`bams[0]` must be a string"));
        assert!(check(&s, &json!({"region": "x", "bams": "one.bam"})).is_ok(), "a lone string is a one-item list");
        assert!(check(&s, &json!({"region": "x", "mode": "c"})).unwrap_err().contains("one of a, b"));
        assert!(check(&s, &json!({"region": "x", "track": 2})).is_ok());
        assert!(check(&s, &json!({"region": "x", "track": true})).is_err());
        assert!(check(&s, &json!({"region": "x", "min_mapq": null})).is_ok(), "null is left out");
        assert!(check(&s, &json!({"region": "x", "other": [1]})).is_ok(), "keys the schema doesn't name are not checked");
    }
}
