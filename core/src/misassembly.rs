use std::{cmp::Ordering, convert::Infallible, str::FromStr};

use serde::Deserialize;

use crate::repeats::Repeat;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Hash)]
pub enum MisassemblyType {
    LowQuality,
    Indel,
    SoftClip,
    Collapse,
    Misjoin,
    FalseDupe,
    RepeatError(Repeat),
    Null,
}

impl MisassemblyType {
    pub fn is_mergeable(&self) -> bool {
        match self {
            MisassemblyType::LowQuality
            | MisassemblyType::SoftClip
            | MisassemblyType::Collapse
            | MisassemblyType::Misjoin => true,
            MisassemblyType::Indel
            | MisassemblyType::FalseDupe
            | MisassemblyType::RepeatError(_)
            | MisassemblyType::Null => false,
        }
    }

    pub fn item_rgb(&self) -> &'static str {
        match self {
            // Purple
            // #800080
            MisassemblyType::Indel => "128,0,128",
            // Teal
            //  #80FFFF
            MisassemblyType::SoftClip => "0,255,255",
            // Pink
            // #FF0080
            MisassemblyType::LowQuality => "255,0,128",
            // Green
            // #00FF00
            MisassemblyType::Collapse => "0,255,0",
            // Orange
            // #FF8000
            MisassemblyType::Misjoin => "255,128,0",
            // Blue
            // #0000FF
            MisassemblyType::FalseDupe => "0,0,255",
            // Scaffold
            // #808080
            MisassemblyType::RepeatError(Repeat::Scaffold) => "128,128,128",
            // Yellow
            // #ECEC00
            MisassemblyType::RepeatError(Repeat::Homopolymer) => "236,236,0",
            // Prussian Blue
            // #003153
            MisassemblyType::RepeatError(Repeat::Dinucleotide) => "0,49,83",
            // Dark green
            // #336600
            MisassemblyType::RepeatError(Repeat::Simple) => "51,102,0",
            // Dirt
            // #9b7653
            MisassemblyType::RepeatError(Repeat::Other) => "155,118,83",
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
            MisassemblyType::RepeatError(Repeat::Scaffold) => "scaffold",
            MisassemblyType::RepeatError(Repeat::Homopolymer) => "homopolymer",
            MisassemblyType::RepeatError(Repeat::Dinucleotide) => "dinucleotide",
            MisassemblyType::RepeatError(Repeat::Simple) => "simple_repeat",
            MisassemblyType::RepeatError(Repeat::Other) => "other_repeat",
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
            "scaffold" => MisassemblyType::RepeatError(Repeat::Scaffold),
            "homopolymer" => MisassemblyType::RepeatError(Repeat::Homopolymer),
            "dinucleotide" => MisassemblyType::RepeatError(Repeat::Dinucleotide),
            "simple_repeat" => MisassemblyType::RepeatError(Repeat::Simple),
            "other_repeat" => MisassemblyType::RepeatError(Repeat::Other),
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
            | (MisassemblyType::Null, MisassemblyType::Null)
            | (MisassemblyType::FalseDupe, MisassemblyType::FalseDupe) => Ordering::Equal,
            // Null/good always less
            (_, MisassemblyType::Null) => Ordering::Greater,
            // Never merge false dupes with others.
            (MisassemblyType::FalseDupe, _) => Ordering::Less,
            (_, MisassemblyType::FalseDupe) => Ordering::Less,
            // Indel and low quality will never replace each other.
            (MisassemblyType::LowQuality, _) => Ordering::Less,
            (MisassemblyType::Indel, _) => Ordering::Less,
            (_, MisassemblyType::Indel) => Ordering::Less,
            // Misjoin should be prioritized over softclip
            (MisassemblyType::SoftClip, MisassemblyType::Misjoin) => Ordering::Less,
            (MisassemblyType::SoftClip, _) => Ordering::Greater,
            // Collapse is greater than misjoin.
            (MisassemblyType::Collapse, MisassemblyType::Misjoin) => Ordering::Greater,
            (MisassemblyType::Collapse, _) => Ordering::Greater,
            // Misjoin always takes priority.
            (MisassemblyType::Misjoin, _) => Ordering::Greater,
            // Never merge repeats
            (_, MisassemblyType::RepeatError(Repeat::Scaffold)) => Ordering::Less,
            (_, MisassemblyType::RepeatError(Repeat::Homopolymer)) => Ordering::Less,
            (_, MisassemblyType::RepeatError(Repeat::Dinucleotide)) => Ordering::Less,
            (_, MisassemblyType::RepeatError(Repeat::Simple)) => Ordering::Less,
            (_, MisassemblyType::RepeatError(Repeat::Other)) => Ordering::Less,
            (MisassemblyType::RepeatError(Repeat::Scaffold), _) => Ordering::Less,
            (MisassemblyType::RepeatError(Repeat::Homopolymer), _) => Ordering::Less,
            (MisassemblyType::RepeatError(Repeat::Dinucleotide), _) => Ordering::Less,
            (MisassemblyType::RepeatError(Repeat::Simple), _) => Ordering::Less,
            (MisassemblyType::RepeatError(Repeat::Other), _) => Ordering::Less,
            // Null/good always less
            (MisassemblyType::Null, _) => Ordering::Less,
        }
    }
}
