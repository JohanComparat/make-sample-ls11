from ls11samples import cli


def test_my_part_covers_everything():
    items = list(range(23))
    parts = [cli.my_part(items, i, 5) for i in range(5)]
    assert sorted(sum(parts, [])) == items
