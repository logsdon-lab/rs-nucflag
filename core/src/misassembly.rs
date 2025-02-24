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
            MisassemblyType::CollapseOther => "yellow",
            MisassemblyType::CollapseVar => "blue",
            MisassemblyType::Collapse => "green",
            MisassemblyType::Misjoin => "orange",
            MisassemblyType::FalseDupe => "purple",
            MisassemblyType::Null => "null",
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
