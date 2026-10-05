import json, pathlib
D = json.load(open(pathlib.Path(__file__).parent.parent / "src" / "aim13_scoping" / "locked_predictions_check.json"))
def test_counts_add_up():
    c = D["counts"]
    assert c["scored"] + c["hedged_called"] + c["unresolved_called"] == c["called"] == 32
    assert len(D["scored"]) == c["scored"] == 5
def test_score():
    assert sum(x["hit"] for x in D["scored"]) == D["counts"]["hits"] == 2
    assert D["counts"]["class_A"] == 0
def test_hit_logic():
    for x in D["scored"]:
        assert x["hit"] == ((x["call"] == "TOLERATED") == x["verdict"].startswith("SAFE"))
