use std::{convert::Infallible, str::FromStr};

use plotters::style::{
    full_palette::{ORANGE, PURPLE},
    Color, RGBAColor, BLUE, GREEN, TRANSPARENT, YELLOW,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MisassemblyType {
    CollapseOther,
    CollapseVar,
    Collapse,
    Misjoin,
    FalseDupe,
    Null,
}

impl From<&MisassemblyType> for RGBAColor {
    fn from(value: &MisassemblyType) -> Self {
        match value {
            MisassemblyType::CollapseOther => YELLOW,
            MisassemblyType::CollapseVar => BLUE,
            MisassemblyType::Collapse => GREEN,
            MisassemblyType::Misjoin => ORANGE,
            MisassemblyType::FalseDupe => PURPLE,
            MisassemblyType::Null => {
                return TRANSPARENT;
            }
        }
        .to_rgba()
        .mix(0.5)
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
