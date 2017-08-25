
types = ["DEL", "INS"]
start_bases = ["T", "C"]
sizes = ["1", "2", "3", "4-10", "10+"]
structures = ["MH", "None", "RPT"]
struc_sizes = ["1", "2", "3+"]

def full_context():
    types = ["DEL", "INS"]
    sizes = ["1", "2", "3", "4-10", "10+"]
    inframe = ["InFrame", "OutFrame"]
    structures = ["RPT", "MH"]
    struc_sizes = ["0", "1", "2", "3", "4-10", "10+"]

    for i in types:
        for j in sizes:
            for k in inframe:
                for l in structures:
                    for m in struc_sizes:
                        print "_".join([i, j , k, l, m])
def part_context():
    for i in types:
        for j in start_bases:
            for k in sizes:
                print "_".join([i, j, k, "None", "X"])
                for t in struc_sizes:
                    if int(k[0]) > 1 and int(t[0]) > 1:
                        print "_".join([i, j , k, "MH", t])
                    print "_".join([i, j, k, "RPT", t])

def tiny_context():
    for i in types:
        for j in start_bases:
            for k in sizes:
                if k == "4-10" or k == "10+":
                    print "_".join([i, j, k, "MH"])
                print "_".join([i, j, k, "NoMH"])

if __name__ == "__main__":
    #part_context()
    tiny_context()
