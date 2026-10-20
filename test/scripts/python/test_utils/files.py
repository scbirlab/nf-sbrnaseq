

def find_repo_root(path, target):
    path = path.resolve()
    for parent in (path.parent, *path.parents):
        if (parent / target).is_file():
            return parent / target
    raise FileNotFoundError("Cannot locate nf-sbrnaseq repository root")
