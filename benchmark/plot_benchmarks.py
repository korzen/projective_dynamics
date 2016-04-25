#!/usr/bin/python3

import sys
import re
import matplotlib.pyplot as plt

class Stats:
    def __init__(self):
        self.max_x = 0
        self.min_x = 0
        self.median = 0
        self.median_abs_dev = 0
        self.mean = 0
        self.stddev = 0
        self.num_samples = 0

    def __repr__(self):
        return self.__str__()
    def __str__(self):
        return "max: {}, min: {}, median: {}, median absdev: {}, stddev: {}".format(
                self.max_x, self.min_x, self.median, self.median_abs_dev, self.stddev)

def parse_benchmarks(fname):
    match_cloth = re.compile("Using cloth of size (\d+)x(\d+)")
    match_mesh = re.compile("Loading mesh (\w+)")
    match_positions = re.compile("n_positions[^\d]*(\d+)")
    match_springs = re.compile("n_springs[^\d]*(\d+)")
    re_stats = {
        "max_x": re.compile(".*max: (\d+\.\d+)"),
        "min_x": re.compile(".*min: (\d+\.\d+)"),
        "median": re.compile(".*median: (\d+\.\d+)"),
        "median_abs_dev": re.compile(".*median abs dev: (\d+\.\d+)"),
        "mean": re.compile(".*mean: (\d+\.\d+)"),
        "stddev": re.compile(".*std dev: (\d+\.\d+)"),
        "num_samples": re.compile(".*# of samples: (\d+)")
    }

    # Benches for the various size cloth and mesh benchmarks
    # A benchmark entry will have:
    #   "n_positions" = # of positions
    #   "n_springs"   = # of springs
    #   "local"       = Stats about local solve
    #   "global"      = Stats about global solve
    #   "full"        = Stats about full solve
    benches = {"cloth": {}, "mesh": {}}
    mesh_type = None
    active_name = None
    active_stats = None
    with open(fname, "r") as f:
        for l in f:
            if l == "pd_benchmark end\n":
                active_name = None
                continue
            if active_name:
                pos = match_positions.match(l)
                if pos:
                    benches[mesh_type][active_name]["n_positions"] = int(pos.group(1))
                    continue
                springs = match_springs.match(l)
                if springs:
                    benches[mesh_type][active_name]["n_springs"] = int(springs.group(1))
                    continue
                if l == "Local Solve stats:\n":
                    active_stats = "local"
                    benches[mesh_type][active_name][active_stats] = Stats()
                elif l == "Global Solve stats:\n":
                    active_stats = "global"
                    benches[mesh_type][active_name][active_stats] = Stats()
                elif l == "Full Solve stats:\n":
                    active_stats = "full"
                    benches[mesh_type][active_name][active_stats] = Stats()
            if active_stats:
                for name, expr in re_stats.items():
                    m = expr.match(l)
                    if m:
                        val = None
                        if name == "num_samples":
                            val = int(m.group(1))
                        else:
                            val = float(m.group(1))
                        setattr(benches[mesh_type][active_name][active_stats], name, val)
                        break

            if not active_name:
                c = match_cloth.match(l)
                mesh = match_mesh.match(l)
                if c:
                    mesh_type = "cloth"
                    active_name = "{}x{}".format(c.group(1), c.group(2))
                    benches[mesh_type][active_name] = {}
                elif mesh:
                    mesh_type = "mesh"
                    active_name = mesh.group(1)
                    benches[mesh_type][active_name] = {}
    return benches

benches = parse_benchmarks(sys.argv[1])

plt.figure()
cloths = ["32x32", "64x64", "128x128", "256x256", "512x512"]
y_val = []
y_err = []
for key in cloths:
    y_val.append(benches["cloth"][key]["global"].mean)
    y_err.append(benches["cloth"][key]["global"].stddev)

plt.errorbar(range(len(cloths)), y_val, yerr=y_err)
plt.xticks(range(len(cloths)), cloths)
plt.ylabel("Time (ms)")
plt.xlabel("Cloth Size")
plt.show()

