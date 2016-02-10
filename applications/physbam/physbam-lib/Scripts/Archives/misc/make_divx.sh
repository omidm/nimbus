#!/bin/bash

bitrate=5000
fps=30
convert_from_ppm=0

while getopts "hb:f:ps:" option
do
	case $option in
		h | "?" )
			echo "Usage: [-b <bitrate>] [-f <fps>] [-s x:y] [-p]";
			exit;;
		b ) bitrate=$OPTARG;;
		f ) fps=$OPTARG;;
        s ) scale="-vf scale=$OPTARG";;
        p ) convert_from_ppm=1;;
	esac
done
shift $(($OPTIND - 1))

echo "bitrate = $bitrate, fps = $fps"
if [ $convert_from_ppm -eq 1 ]; then
    echo "converting files from ppm"
    for i in *.ppm; do
        echo "converting $i"
        convert $i `basename $i .ppm`.png;
    done
fi

mencoder mf://*.png -mf fps=$fps -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=$bitrate -ffourcc XVID $scale -o "movie_$bitrate.avi"

if [ $convert_from_ppm -eq 1 ]; then
    rm *.png
fi
