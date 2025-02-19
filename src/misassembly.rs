use plotters::style::{
    full_palette::{ORANGE, PURPLE, TEAL},
    Color, RGBAColor, BLUE, GREEN, YELLOW,
};

#[derive(Debug, Clone)]
pub enum Misassembly {
    CollapseOther,
    CollapseVar,
    Collapse,
    Misjoin,
    FalseDupe,
    Het,
}

impl From<Misassembly> for RGBAColor {
    fn from(value: Misassembly) -> Self {
        match value {
            Misassembly::CollapseOther => YELLOW,
            Misassembly::CollapseVar => BLUE,
            Misassembly::Collapse => GREEN,
            Misassembly::Misjoin => ORANGE,
            Misassembly::FalseDupe => PURPLE,
            Misassembly::Het => TEAL,
        }
        .to_rgba()
        .mix(0.5)
    }
}
