use std::{convert::Infallible, str::FromStr};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MisassemblyType {
    CollapseOther,
    CollapseVar,
    Collapse,
    Misjoin,
    FalseDupe,
    Null,
}

impl MisassemblyType {
    pub fn item_rgb(&self) -> &'static str {
        match self {
            MisassemblyType::CollapseOther => "255,255,0",
            MisassemblyType::CollapseVar => "0,0,255",
            MisassemblyType::Collapse => "0,255,0",
            MisassemblyType::Misjoin => "255,165,0",
            MisassemblyType::FalseDupe => "128,0,128",
            MisassemblyType::Null => "0,0,0",
        }
    }
}

impl FromStr for MisassemblyType {
    type Err = Infallible;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "collapse_other" => MisassemblyType::CollapseOther,
            "collapse_var" => MisassemblyType::CollapseVar,
            "misjoin" => MisassemblyType::Misjoin,
            "collapse" => MisassemblyType::Collapse,
            "false_dupe" => MisassemblyType::FalseDupe,
            _ => MisassemblyType::Null,
        })
    }
}
