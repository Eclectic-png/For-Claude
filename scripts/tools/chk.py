import bpy, bmesh
for nm in ("Hips","AN_AnalCanal","AN_Rectum_LowerAmpulla","AN_Rectum"):
    o=bpy.data.objects.get(nm)
    if not o: print("CHK",nm,"missing"); continue
    bm=bmesh.new(); bm.from_mesh(o.data)
    print("CHK",nm,len(bm.verts),"boundary",sum(e.is_boundary for e in bm.edges),"nonmanifold",sum(not e.is_manifold and not e.is_boundary for e in bm.edges))
print("CHK engine",bpy.context.scene.render.engine, [o.name for o in bpy.data.objects if o.type in('CAMERA','LIGHT')])
