from video import *
import sys

header_color='#0099cc'
slide("title.png",48,[(50,"A Novel Algorithm for",header_color),
                      (50,"Incompressible Flow",header_color),
                      (50,"Using Only a Coarse",header_color),
                      (50,"Grid Projection",header_color),
                      " ",
                      (36,"paperid 0046"),
                      (36,"SIGGRAPH 2010")])


slide("static_sphere.png",36,
      [(54,"Static Sphere",header_color),
       " ",
       "Smoke with a static sphere"])

slide("static_sphere_base.png",36,
      [(54,"Static Sphere",header_color),
       " ",
       "Base Simulation",
       "Resolution: 64x128x64"])

slide("static_sphere_coarse1.png",36,
      [(54,"Static Sphere",header_color),
       " ",
       "Fine: 64x128x64",
       "Coarse: 32x64x32"])

slide("static_sphere_coarse2.png",36,
      [(54,"Static Sphere",header_color),
       " ",
       "Fine: 64x128x64",
       "Coarse: 16x32x16",
       " ",
       "(with Kolmolgorov noise)"])

slide("static_sphere_refine1.png",36,
      [(54,"Static Sphere",header_color),
       " ",
       "Fine: 128x256x128",
       "Coarse: 64x128x64"])

slide("static_sphere_refine2.png",36,
      [(54,"Static Sphere",header_color),
       " ",
       "Fine: 256x512x256",
       "Coarse: 64x128x64"])

slide("static_sphere_refine3.png",36,
      [(54,"Static Sphere",header_color),
       " ",
       "Fine: 512x1024x512",
       "Coarse: 64x128x64"])

slide("static_sphere_compare.png",36,
      [(54,"Static Sphere",header_color),
       " ",
       "Side by side comparison",
       " ",
       "3 refinement levels"])

slide("motion_sphere.png",36,
      [(54,"Moving Sphere",header_color),
       " ",
       "Smoke with a moving sphere"])

slide("motion_sphere_base.png",36,
      [(54,"Moving Sphere",header_color),
       " ",
       "Base Simulation",
       "Resolution: 64x128x64"])

slide("motion_sphere_coarse1.png",36,
      [(54,"Moving Sphere",header_color),
       " ",
       "Fine: 64x128x64",
       "Coarse: 32x64x32"])

slide("motion_sphere_coarse2.png",36,
      [(54,"Moving Sphere",header_color),
       " ",
       "Fine: 64x128x64",
       "Coarse: 16x32x16",
       " ",
       "(with Kolmolgorov noise)",
       ])

slide("motion_sphere_refine1.png",36,
      [(54,"Moving Sphere",header_color),
       " ",
       "Fine: 128x256x128",
       "Coarse: 64x128x64"])

slide("motion_sphere_refine2.png",36,
      [(54,"Moving Sphere",header_color),
       " ",
       "Fine: 256x512x256",
       "Coarse: 64x128x64"])

slide("motion_sphere_refine3.png",36,
      [(54,"Moving Sphere",header_color),
       " ",
       "Fine: 512x1024x512",
       "Coarse: 64x128x64"])

slide("motion_sphere_compare.png",36,
      [(54,"Moving Sphere",header_color),
       " ",
       "Side by side comparison",
       " ",
       "3 refinement levels"])

slide("water.png",36,
      [(54,"Water",header_color),
       " ",
       "Water falling into a box"])

slide("water_base.png",36,
      [(54,"Water",header_color),
       " ",
       "Base Simulation",
       "Resolution: 64x128x64"])

slide("water_coarse1.png",36,
      [(54,"Water",header_color),
       " ",
       "Fine: 128x128x128",
       "Coarse: 64x64x64"])

slide("water_coarse2.png",36,
      [(54,"Water",header_color),
       " ",
       "Fine: 128x128x128",
       "Coarse: 32x32x32"])

slide("water_refine1.png",36,
      [(54,"Water",header_color),
       " ",
       "Fine: 256x256x256",
       "Coarse: 128x128x128"])

slide("water_refine2.png",36,
      [(54,"Water",header_color),
       " ",
       "Fine: 512x512x512",
       "Coarse: 128x128x128"])

slide("water_compare.png",36,
      [(54,"Water",header_color),
       " ",
       "Side by side comparison",
       " ",
       "2 refinement levels"])

slide("compare.png",36,
      ["base                   our",
      "simulation           method"])

slide("overview.png",36,
      ["Algorithmic Overview"])

final_render_base="/solver/vol0/refinement7/render"

video_dir="refinement_video_tmp"
if not os.path.exists(video_dir):
    os.mkdir(video_dir)
#if os.path.isdir(video_dir):
#    print "%s already exists... delete? (y/n) "%video_dir,
#    c=sys.stdin.readline()
#    if c=="y\n": shutil.rmtree(video_dir)
#    else: sys.exit(1)
shutil.rmtree(video_dir)    
video=VIDEO(video_dir)

using_lengths=True
testing_titles=False
skip_interval=1
hi_res=1

def caption_length(time=3.5):
    if not using_lengths and testing_titles:
        return 1
    else:
        return int(video.fps*time)

video.add_frame("title.png",caption_length(7./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_motion/refine_triple_large_pngs",1,70,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_motion/refine_triple_pngs",1,70,step=skip_interval)

image_dir="/home/zhw/public_html/Refinement_Instruction_Video/output/"
frame_rate=24.
if hi_res:
    video.add_frame(final_render_base+"/3D_Smoke_motion/refine_triple_large_pngs/frame.00070.png",int(.5*frame_rate))
    video.add_blending(final_render_base+"/3D_Smoke_motion/refine_triple_large_pngs/frame.00070.png","overview.png",int(2.*frame_rate))
else:
    video.add_frame(final_render_base+"/3D_Smoke_motion/refine_triple_pngs/frame.00070.png",int(.5*frame_rate))
    video.add_blending(final_render_base+"/3D_Smoke_motion/refine_triple_pngs/frame.00070.png","overview.png",int(2.*frame_rate))
video.add_frame("overview.png",caption_length(3./skip_interval))
video.add_frame(image_dir+"img0.png",int(2.*frame_rate))
video.add_blending(image_dir+"img0.png",image_dir+"img1.png",int(2.*frame_rate))
video.add_frame(image_dir+"img1.png",int(2.*frame_rate))
video.add_blending(image_dir+"img1.png",image_dir+"img2.png",int(2.*frame_rate))
video.add_frame(image_dir+"img2.png",int(2.*frame_rate))
video.add_blending(image_dir+"img2.png",image_dir+"img3.png",int(2.*frame_rate))
video.add_frame(image_dir+"img3.png",int(2.*frame_rate))
video.add_blending(image_dir+"img3.png",image_dir+"img4.png",int(2.*frame_rate))
video.add_frame(image_dir+"img4.png",int(2.*frame_rate))
video.add_blending(image_dir+"img4.png",image_dir+"img5.png",int(2.*frame_rate))
video.add_frame(image_dir+"img5.png",int(2.*frame_rate))

video.add_frame("static_sphere.png",caption_length(5./skip_interval))
video.add_frame("static_sphere_base.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_static/base_64_large_pngs",1,120,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_static/base_64_pngs",1,120,step=skip_interval)

video.add_frame("static_sphere_coarse1.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_static/coarsen_once_large_pngs",1,120,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_static/coarsen_once_pngs",1,120,step=skip_interval)

video.add_frame("static_sphere_coarse2.png",caption_length(5./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_static/coarsen_twice_large_pngs",1,120,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_static/coarsen_twice_pngs",1,120,step=skip_interval)

video.add_frame("static_sphere_refine1.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_static/refine_once_large_pngs",1,120,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_static/refine_once_pngs",1,120,step=skip_interval)

video.add_frame("static_sphere_refine2.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_static/refine_twice_large_pngs",1,120,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_static/refine_twice_pngs",1,120,step=skip_interval)

video.add_frame("static_sphere_refine3.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_static/refine_triple_large_pngs",1,120,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_static/refine_triple_pngs",1,120,step=skip_interval)

video.add_frame("static_sphere_compare.png",caption_length(5./skip_interval))
video.add_frame("compare.png",caption_length(2./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_static/compare_large_pngs",1,120,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_static/compare_small_pngs",1,120,step=skip_interval)

video.add_frame("motion_sphere.png",caption_length(5./skip_interval))
video.add_frame("motion_sphere_base.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_motion/base_64_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_motion/base_64_pngs",1,200,step=skip_interval)

video.add_frame("motion_sphere_coarse1.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_motion/coarsen_once_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_motion/coarsen_once_pngs",1,200,step=skip_interval)

video.add_frame("motion_sphere_coarse2.png",caption_length(5./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_motion/coarsen_twice_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_motion/coarsen_twice_pngs",1,200,step=skip_interval)

video.add_frame("motion_sphere_refine1.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_motion/refine_once_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_motion/refine_once_pngs",1,200,step=skip_interval)

video.add_frame("motion_sphere_refine2.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_motion/refine_twice_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_motion/refine_twice_pngs",1,200,step=skip_interval)

video.add_frame("motion_sphere_refine3.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_motion/refine_triple_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_motion/refine_triple_pngs",1,200,step=skip_interval)

video.add_frame("motion_sphere_compare.png",caption_length(5./skip_interval))
video.add_frame("compare.png",caption_length(2./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Smoke_motion/compare_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Smoke_motion/compare_small_pngs",1,200,step=skip_interval)

video.add_frame("water.png",caption_length(5./skip_interval))
video.add_frame("water_base.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Water/128_base_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Water/128_base_pngs",1,200,step=skip_interval)

video.add_frame("water_coarse1.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Water/128_coarse_1_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Water/128_coarse_1_pngs",1,200,step=skip_interval)

video.add_frame("water_coarse2.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Water/128_coarse_2_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Water/128_coarse_2_pngs",1,200,step=skip_interval)

video.add_frame("water_refine1.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Water/128_refined_1_surface_large_pngs",1,200,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Water/128_refined_1_surface_pngs",1,200,step=skip_interval)

video.add_frame("water_refine2.png",caption_length(4./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Water/128_refined_2_large_reintialized_pngs",1,162,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Water/128_refined_2_pngs",1,162,step=skip_interval)

video.add_frame("water_compare.png",caption_length(5./skip_interval))
video.add_frame("compare.png",caption_length(2./skip_interval))
if not testing_titles:
    if hi_res: video.add_directory(final_render_base+"/3D_Water/compare_large_pngs",1,162,step=skip_interval)
    else: video.add_directory(final_render_base+"/3D_Water/compare_small_pngs",1,162,step=skip_interval)

#video.add_frame("end.png",caption_length(2./skip_interval))

video.make_movie('refinement')


