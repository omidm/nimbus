#!/usr/bin/python

from video import *
import Image
import sys

header_color='#ffaa00'

def Make_Slides():
    slide("title.png",40,
            [(40,"Practical Animation of",header_color), 
             (40,"Compressible Flow",header_color), 
             (40,"for Shock Waves",header_color), 
             (40,"and",header_color), 
             (40,"Related Phenomena",header_color), 
             "",
             (30,"paperid 1030"),
             (30,"SCA 2010")])

    slide("shock_hit_heavy_wall.png",36,
            [(45,"Spherical Shock",header_color),
             (45,"Hitting a Heavy Solid",header_color),
             "",
             "Notice Strong Shock Reflection",
             "Off the Wall"
             ])

    slide("shock_hit_light_wall.png",36,
            [(45,"Spherical Shock",header_color),
             (45,"Hitting a Light Solid",header_color),
             "",
             "Shock Weakly Reflects",
             "Notice Secondary Shocks",
             "Generated from Collision"
             ])

    slide("obrien.png",36,
            [(45,"Pressure Profile of",header_color),
             (45,"Shock Hitting",header_color),
             (45,"a Static Solid",header_color),
             "",
             "2-D Slice of 3-D Simulation"
             "",
             "Comparison with Yngve et al.",
             "(Animating Explosions",
             "SIGGRAPH 2000)"
             ])

    slide("trinity.png",36,
            [(45,"Trinity Test",header_color),
             "",
             ("Initial Pressure=9.4e10 Pa"),
             ("Initial Temperature=2.6e8 K"),
             ])

    slide("shock_hit_smoke1.png",36,
            [(45,"Smoke Plume",header_color),
             (45,"Hit by a Shock",header_color),
             ])

    slide("shock_hit_smoke2.png",36,
            [(45,"Shock Introduced"),
             (45,"at Left"),
             "",
             "(Slow Motion at 25000X)",
             "Smoke Barely Moves Vertically",
             "at this Time Scale"
             ])

    slide("stack1.png",36,
            [(45,"Stack of Rigid Bodies",header_color),
             (45,"Hit by a Shock",header_color),
             ])

    slide("stack2.png",36,
            ["Slow Motion (50X)",
             "to Show Shock Front "
             ])
   
    # try compositing
    slide("stack3.png",36,
            ["A Little Faster...",
             ])

    slide("cannon1.png",36,
            [(45,"Cannon Fired At",header_color),
             (45,"Bunny",header_color),
             ])

    slide("cannon2.png",36,
            ["Slow Motion (4X)",
             "Notice Secondary Shock",
             "Created by",
             "Supersonic Cannon Ball"
             ])

    slide("deformable1.png",36,
            [(45,"Deformable Ball",header_color),
             (45,"Hit By",header_color),
             (45,"Spherical Shock",header_color),
             ])

    slide("deformable2.png",36,
            ["Different View",
             "Without Rendering",
             "Soot and Blackbody",
             "to Show Ball Deformation"
             ])

    slide("spherical_four_walls1.png",36,
            [(45,"Shock in a",header_color),
             (45,"Detonation Chamber",header_color),
             "",
             "Smoothly Transitions",
             "to Incompressible Halfway",
             ])

    slide("spherical_four_walls2.png",36,
             ["Slow Motion",
             "To Show Shocks",
             "",
             "(Annotation Shows",
             "State of Transition)"
             ])

    slide("spherical_four_walls3.png",36,
             ["A Little Faster...",
             "",
             "(Annotation Shows",
             "State of Transition)"
             ])

    slide("spherical_four_walls_fracture1.png",36,
            [(45,"Shock in a",header_color),
             (45,"Detonation Chamber",header_color),
             (45,"With a Fragile Wall",header_color),
             "",
             "Smoothly Transitions",
             "to Incompressible Halfway",
             ])

    slide("spherical_four_walls_fracture2.png",36,
             ["Slow Motion",
             "To Show Shocks",
             ])

    slide("spherical_four_walls_fracture3.png",36,
             ["A Little Faster...",
             ])

    slide("spherical_four_walls_fracture_annotate1.png",36,
             ["Slow Motion",
             "To Show Shocks",
             "",
             "(Annotation Shows",
             "State of Transition)"
             ])

    slide("spherical_four_walls_fracture_annotate2.png",36,
             ["A Little Faster...",
             "",
             "(Annotation Shows",
             "State of Transition)"
             ])

    slide("spherical_four_walls_fracture_zoom.png",36,
             ["Close Up View of Fracture"
             ])

    slide("spherical_four_walls_fracture_zoom_nofire.png",36,
             ["Close Up View of Fracture",
             "Without Rendering",
             "Soot and Blackbody",
             "to Show Fracture"
             ])

def Make_Video(video_dir):
    video=VIDEO(video_dir)
    title_length=int(video.fps*10)
    caption_length=int(video.fps*4)
    one_second=video.fps
    final_render_base="/n/self/disk2/solid_compressible_coupling/FINAL_RENDERS"

    video.add_frame("title.png",caption_length)
   
    # Heavy wall
    video.add_frame("shock_hit_heavy_wall.png",one_second*6)
    video.add_directory(final_render_base+"/shock_hit_heavy_wall/renders",0,2000,4)

    # Light wall
    video.add_frame("shock_hit_light_wall.png",int(one_second*7.5))
    video.add_directory(final_render_base+"/shock_hit_light_wall/renders",0,2000,4)

    # Obrien 
    video.add_frame("obrien.png",one_second*5)
    video.add_directory(final_render_base+"/obrien/renders",0,264,1)

    # Trinity 
    video.add_frame("trinity.png",one_second*5)
    video.add_directory(final_render_base+"/trinity/renders_black_body",0,104,1)

    # Shock hit smoke 
    video.add_frame("shock_hit_smoke1.png",one_second*3)
    video.add_directory(final_render_base+"/shock_hit_smoke/renders_smoke_only",0,50,1)

    video.composite_slide_on_current_frame("shock_hit_smoke2.png","shock_hit_smoke2_composite.png",one_second*7)
    video.add_directory(final_render_base+"/shock_hit_smoke/renders_shock_hitting",0,350,1)

    # stack
    video.add_frame("stack1.png",one_second*4)
    video.add_directory(final_render_base+"/stack/renders",1,3991,50)

    video.add_frame("stack2.png",caption_length)
    video.add_directory(final_render_base+"/stack/renders",1,601,1)

    video.composite_slide_on_current_frame("stack3.png","stack3_composite.png",int(one_second*2.5))
    video.add_directory(final_render_base+"/stack/renders",601,4000,10)

    # Cannon
    video.add_frame("cannon1.png",one_second*4)
    video.add_directory(final_render_base+"/cannon_bunny/renders",0,700,8)

    video.add_frame("cannon2.png",int(one_second*5.5))
    video.add_directory(final_render_base+"/cannon_bunny/renders",0,700,2)

    # Deformable
    video.add_frame("deformable1.png",one_second*5)
    video.add_directory(final_render_base+"/deformable/renders",0,2990,10)

    video.add_frame("deformable2.png",one_second*5)
    video.add_directory(final_render_base+"/deformable_shock_only/renders",0,2990,10)

    # Four walls
    video.add_frame("spherical_four_walls1.png",caption_length)
    video.add_directory(final_render_base+"/spherical_4_walls_highres/renders_bottomview",0,3500,10)

    video.add_frame("spherical_four_walls2.png",caption_length)
    video.add_directory(final_render_base+"/spherical_4_walls_highres/renders_bottomview_composite",0,200,1)

    video.composite_slide_on_current_frame("spherical_four_walls3.png","spherical_four_walls3_composite.png",one_second*3)
    video.add_directory(final_render_base+"/spherical_4_walls_highres/renders_bottomview_composite",200,3500,10)

    # Fracture
    video.add_frame("spherical_four_walls_fracture1.png",caption_length)
    video.add_directory(final_render_base+"/spherical_4_walls_fracture_highres/renders_bottomview",0,2830,10)

    video.add_frame("spherical_four_walls_fracture_annotate1.png",caption_length)
    video.add_directory(final_render_base+"/spherical_4_walls_fracture_highres/renders_bottomview_composite",0,200,1)

    video.composite_slide_on_current_frame("spherical_four_walls_fracture_annotate2.png","spherical_four_walls_fracture3_composite.png",one_second*3)
    video.add_directory(final_render_base+"/spherical_4_walls_fracture_highres/renders_bottomview_composite",200,2830,10)

    video.add_frame("spherical_four_walls_fracture_zoom_nofire.png",one_second*5)
    video.add_directory(final_render_base+"/spherical_4_walls_fracture_highres/renders_zoomed",0,120,1)

    # Make the video
    video.make_movie('compressible')

def main():
    video_dir="compressible_video"
    if not os.path.exists(video_dir):
        os.mkdir(video_dir)

    Make_Slides()
    Make_Video(video_dir)

if __name__=="__main__":
    main()

