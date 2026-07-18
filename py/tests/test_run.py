import tempfile
import tomlkit
import py_nucflag

BAM = "../core/test/standard/input/aln_1.bam"
BED = "../core/test/standard/input/aln_1.bed"
CFG = "tests/cfg.toml"
ITV = (3021508, 3021608, "haplotype2-0000133")

def test_run_nucflag():
    res = py_nucflag.run_nucflag(BAM, bed=BED)[0]
    assert res.pileup.shape == (5_669_966, 13)
    assert res.regions.shape == (58, 11) 

def test_run_nucflag_itv():
    res = py_nucflag.run_nucflag_itv(BAM, itv=ITV)
    assert res.pileup.shape == (101, 13)
    assert res.regions.shape == (1, 11)

def test_get_regions_no_bed():
    itvs = py_nucflag.get_regions(BAM)
    assert len(itvs) == 1263

def test_get_regions_w_bed():
    itvs = py_nucflag.get_regions(BAM, bed=BED)
    assert itvs == [(3021508, 8691473, 'haplotype2-0000133')]

def test_get_regions_w_invalid_bed():
    with tempfile.NamedTemporaryFile("wt") as fh:
        print("chr?\t0\t100", file=fh)
        fh.flush()
        res = py_nucflag.get_regions(BAM, bed=fh.name)
        assert len(res) == 0

def test_get_config_from_preset():
    res = tomlkit.parse(py_nucflag.get_config_from_preset("hifi"))
    res["repeat"]["check_types"].sort()
    with open(CFG, "rt") as fh:
        exp = tomlkit.parse(fh.read())
        exp["repeat"]["check_types"].sort()
    assert res == exp
