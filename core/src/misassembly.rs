use std::{cmp::Ordering, convert::Infallible, str::FromStr};

use serde::Deserialize;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Hash)]
pub enum MisassemblyType {
    LowQuality,
    Indel,
    SoftClip,
    Collapse,
    Misjoin,
    FalseDupe,
    Null,
}

impl MisassemblyType {
    pub fn item_rgb(&self) -> &'static str {
        match self {
            // Purple
            MisassemblyType::Indel => "128,0,128",
            // Teal
            MisassemblyType::SoftClip => "0,255,255",
            // Pink
            MisassemblyType::LowQuality => "255,0,128",
            // Green
            MisassemblyType::Collapse => "0,255,0",
            // Orange
            MisassemblyType::Misjoin => "255,165,0",
            // Blue
            MisassemblyType::FalseDupe => "0,0,255",
            MisassemblyType::Null => "0,0,0",
        }
    }
}

impl From<MisassemblyType> for &'static str {
    fn from(value: MisassemblyType) -> Self {
        match value {
            MisassemblyType::LowQuality => "low_quality",
            MisassemblyType::Indel => "indel",
            MisassemblyType::SoftClip => "softclip",
            MisassemblyType::Collapse => "collapse",
            MisassemblyType::Misjoin => "misjoin",
            MisassemblyType::FalseDupe => "false_dupe",
            MisassemblyType::Null => "null",
        }
    }
}

impl FromStr for MisassemblyType {
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "low_quality" => MisassemblyType::LowQuality,
            "indel" => MisassemblyType::Indel,
            "softclip" => MisassemblyType::SoftClip,
            "misjoin" => MisassemblyType::Misjoin,
            "collapse" => MisassemblyType::Collapse,
            "false_dupe" => MisassemblyType::FalseDupe,
            _ => MisassemblyType::Null,
        })
    }
}

impl PartialOrd for MisassemblyType {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MisassemblyType {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            // Equal if same.
            (MisassemblyType::LowQuality, MisassemblyType::LowQuality)
            | (MisassemblyType::Indel, MisassemblyType::Indel)
            | (MisassemblyType::SoftClip, MisassemblyType::SoftClip)
            | (MisassemblyType::Collapse, MisassemblyType::Collapse)
            | (MisassemblyType::Misjoin, MisassemblyType::Misjoin)
            | (MisassemblyType::FalseDupe, MisassemblyType::FalseDupe) => Ordering::Equal,
            // Indel and low quality will never replace each other.
            (MisassemblyType::LowQuality, _) => Ordering::Less,
            (MisassemblyType::Indel, _) => Ordering::Less,
            // Misjoin should be prioritized over softclip
            (MisassemblyType::SoftClip, MisassemblyType::Misjoin) => Ordering::Less,
            (MisassemblyType::SoftClip, _) => Ordering::Greater,
            // Collapse is less than misjoin.
            (MisassemblyType::Collapse, MisassemblyType::Misjoin) => Ordering::Less,
            (MisassemblyType::Collapse, _) => Ordering::Greater,
            // Misjoin and false dupe always takes priority.
            (MisassemblyType::Misjoin, _) => Ordering::Greater,
            (MisassemblyType::FalseDupe, _) => Ordering::Greater,
            (MisassemblyType::Null, _) => unreachable!("Null misassembly type."),
        }
    }
}
