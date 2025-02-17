#[derive(Debug, Clone)]
pub enum Misassembly {
    CollapseOther,
    CollapseVar,
    Collapse,
    Misjoin,
    FalseDupe,
    Het,
}
