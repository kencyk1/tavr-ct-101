#!/bin/bash
cd /home/claude/seg
for part in organs vertebrae cardiac muscles ribs; do
  names=$(python3 -c "from totalsegmentator.map_to_binary import class_map_5_parts as m;print(' '.join(m['class_map_part_$part'].values()))")
  TotalSegmentator -i ct1349.nii.gz -o part_$part.nii --ml -nr 1 -ns 1 -fs --roi_subset $names > part_$part.log 2>&1
  echo "$part exit $?"
  pkill -9 -f "python3.13 -c from"; sleep 2
done
