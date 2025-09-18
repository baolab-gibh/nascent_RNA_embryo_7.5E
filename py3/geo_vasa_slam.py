def convert_pos(adt, into="cell_type"):
    '''Convert sample name into cell type, regions, or layers.'''
    _pos = ['A', 'P', 'EA', 'EP', 'L', 'R', 'MA', 'MP', 'O']
    if into == "regions":
        return adt.obs.index.str.extract("[0-9]+([A-Z]{1,2})").loc[:, 0].tolist()
    elif into == "layers":
        return adt.obs.index.str.extract("(NC_[0-9]+|[0-9]+)[A-Z]{1,2}").loc[:, 0].tolist()
    else:
        if into == "pseudotime":
            #_target = [3, 3, 1, 1, 3, 3, 2, 2, 0]
            _target = [1, 1, 1, 1, 1, 1, 1, 1, 0]
        else:
            _target = ["Ectoderm", "Ectoderm", "Endoderm", "Endoderm", "Ectoderm", "Ectoderm", "Mesoderm", "Mesoderm", "Other"]

        pseudo_time = dict(zip(_pos, _target))
        regions = [x[0] if x else "O" for x in adt.obs.index.str.findall("[0-9]+([A-Z]{1,2})")]
        return [pseudo_time[x] for x in regions]


