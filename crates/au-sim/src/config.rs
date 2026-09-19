//! Configuration — because constants are data, not code.
//!
//! VISION, Programming Philosophy: "Data-driven systems over hardcoded
//! behavior." Physical constants are the first place that principle gets
//! tested. If gravity is a literal in a `.rs` file, then "what if the constants
//! were different?" — a question this project is *supposed* to be able to
//! ask — requires a recompile.
//!
//! Tiny hand-written parser rather than a TOML crate, for the same reason
//! au-core has no dependencies: this file's contents change the outcome of the
//! universe, and I want to know exactly how they are parsed.
//!
//! Format: `key = value`, `#` comments, blank lines ignored.

use std::collections::BTreeMap;
use std::fmt;

#[derive(Debug)]
pub enum ConfigError {
    Syntax { line: usize, text: String },
    Missing(String),
    BadValue { key: String, value: String, want: &'static str },
}

impl fmt::Display for ConfigError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ConfigError::Syntax { line, text } => {
                write!(f, "config line {}: expected `key = value`, got `{}`", line, text)
            }
            ConfigError::Missing(k) => write!(f, "config: missing required key `{}`", k),
            ConfigError::BadValue { key, value, want } => {
                write!(f, "config: `{}` = `{}` is not a valid {}", key, value, want)
            }
        }
    }
}

impl std::error::Error for ConfigError {}

#[derive(Clone, Debug, Default)]
pub struct Config {
    values: BTreeMap<String, String>,
}

impl Config {
    pub fn parse(src: &str) -> Result<Config, ConfigError> {
        let mut values = BTreeMap::new();
        for (i, raw) in src.lines().enumerate() {
            let line = raw.split('#').next().unwrap_or("").trim();
            if line.is_empty() {
                continue;
            }
            let (k, v) = line.split_once('=').ok_or_else(|| ConfigError::Syntax {
                line: i + 1,
                text: line.to_string(),
            })?;
            values.insert(k.trim().to_string(), v.trim().to_string());
        }
        Ok(Config { values })
    }

    fn raw(&self, key: &str) -> Result<&str, ConfigError> {
        self.values
            .get(key)
            .map(|s| s.as_str())
            .ok_or_else(|| ConfigError::Missing(key.to_string()))
    }

    pub fn u64(&self, key: &str) -> Result<u64, ConfigError> {
        let v = self.raw(key)?;
        v.parse().map_err(|_| ConfigError::BadValue {
            key: key.into(),
            value: v.into(),
            want: "integer",
        })
    }

    pub fn f64(&self, key: &str) -> Result<f64, ConfigError> {
        let v = self.raw(key)?;
        v.parse().map_err(|_| ConfigError::BadValue {
            key: key.into(),
            value: v.into(),
            want: "number",
        })
    }

    pub fn bool(&self, key: &str) -> Result<bool, ConfigError> {
        let v = self.raw(key)?;
        match v {
            "true" | "yes" | "1" => Ok(true),
            "false" | "no" | "0" => Ok(false),
            _ => Err(ConfigError::BadValue { key: key.into(), value: v.into(), want: "bool" }),
        }
    }

    pub fn u64_or(&self, key: &str, default: u64) -> u64 {
        self.u64(key).unwrap_or(default)
    }
    pub fn f64_or(&self, key: &str, default: f64) -> f64 {
        self.f64(key).unwrap_or(default)
    }
    pub fn bool_or(&self, key: &str, default: bool) -> bool {
        self.bool(key).unwrap_or(default)
    }
    /// Absent and present-but-unset are different things. A missing
    /// `bottom_temp_k` means "no fixed-temperature wall here", not "zero kelvin".
    pub fn f64_opt(&self, key: &str) -> Option<f64> {
        self.f64(key).ok()
    }
    pub fn str_or<'a>(&'a self, key: &str, default: &'a str) -> &'a str {
        self.raw(key).unwrap_or(default)
    }

    pub fn set(&mut self, key: &str, value: &str) {
        self.values.insert(key.into(), value.into());
    }

    /// Overlay another config file on top of this one. Later wins.
    ///
    /// This is how the validation fixtures (`data/reference_materials.kv`) get
    /// in without the engine ever shipping a substance of its own.
    pub fn merge(&mut self, src: &str) -> Result<(), ConfigError> {
        let other = Config::parse(src)?;
        for (k, v) in other.values {
            self.values.insert(k, v);
        }
        Ok(())
    }

    pub fn as_map(&self) -> &BTreeMap<String, String> {
        &self.values
    }

    pub fn iter(&self) -> impl Iterator<Item = (&String, &String)> {
        self.values.iter()
    }
}

/// The defaults, embedded so the engine runs with no files present.
/// `data/constants.kv` overrides them.
pub const DEFAULT_CONFIG: &str = include_str!("../../../data/constants.kv");
