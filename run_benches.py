#!/usr/bin/python3

import os

start_cloth = 32
end_cloth = 512

cloth = start_cloth
out_name = "viennacl_jacobi"
do_cloth = True
if do_cloth:
    while cloth <= end_cloth:
        os.system("./pd_benchmark 1000 240 --size {} {} | tee -a benchmark/{}.txt".format(cloth, cloth, out_name))
        cloth = cloth * 2

meshes = ["chihuahua", "ant", "rhino"]
do_meshes = False
if do_meshes:
    for m in meshes:
        os.system("./pd_benchmark 1000 240 --mesh ./mesh/{}.bmesh | tee -a benchmark/{}.txt".format(m, out_name))

