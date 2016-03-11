bl_info = {
        "name":     "Projective Dynamics Export",
        "author":   "Pavol Klacansky and Will Usher",
        "version":  (0, 1),
        "blender":  (2, 76, 0),
        "location": "File > Export > Projective Dynamics (.json)",
        "category": "Import-Export",
}


import bpy
from bpy_extras.io_utils import ExportHelper
import itertools
import json


def convert_matrix(mat):
        return list(itertools.chain.from_iterable([row[:] for row in mat]))


# TODO: assume mesh was triangulated
def convert_mesh(mesh):
        vertices    = [v.co[:] for v in mesh.data.vertices]
        indices     = [s.vertices[:] for s in list(mesh.data.polygons)]
        attachments = []

        cloth_modifiers = [m for m in mesh.modifiers if m.type == "CLOTH"]
        if cloth_modifiers:
            cloth = cloth_modifiers[0]
            group = cloth.settings.vertex_group_mass
            if cloth.settings.use_pin_cloth and group != "":
                index = mesh.vertex_groups[group].index
                attachments = [{ "i": i, "position": vertices[i] }  for i, v in enumerate(mesh.data.vertices) for g in v.groups if g.group == index]

        # in case of triangle mesh, springs are just edge set
        springs = [{ "i": list(e.vertices), "rest_length": (mesh.data.vertices[e.vertices[0]].co - mesh.data.vertices[e.vertices[1]].co).length } for e in mesh.data.edges]

        return {
                "vertices":    vertices,
                "indices":     indices,
                "attachments": attachments,
                "springs":     springs
        }


def convert_object(obj):
        data = {
#                "matrix_world": convert_matrix(obj.matrix_world),
                "name": obj.name,
        }

        # we support only meshes
        if obj.type != "MESH":
                return

        data["type"] = "mesh"
        data["data"] = convert_mesh(obj)

        return data


def export(context, props, filepath):
        scene = context.scene

        exp_objects = [convert_object(obj) for obj in scene.objects]
        exp_objects = list(filter(None.__ne__, exp_objects))

        exp_scene = {
                "objects": exp_objects
        }

        with open(filepath, "w") as fp:
                json.dump(exp_scene, fp, indent=8, sort_keys=True)


class PDExport(bpy.types.Operator, ExportHelper):
        bl_idname = "export_scene.pd_export"
        bl_label = "Export PD"
        filename_ext = ".json"

        def execute(self, context):
                export(context, self.properties, self.filepath)
                return {"FINISHED"}


def create_menu(self, context):
        self.layout.operator(PDExport.bl_idname,
                             text="Projective Dynamics (.json)")


def register():
        bpy.utils.register_module(__name__)
        bpy.types.INFO_MT_file_export.append(create_menu)


def unregister():
        bpy.utils.unregister_module(__name__)
        bpy.types.INFO_MT_file_export.remove(create_menu)


if __name__ == "__main__":
        register()
