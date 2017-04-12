# this script is meant to be run from a machine on UW campus

output_dir=/astro/store/pogo4/danielsf/mlt_flares_shorter/

nohup nice python assign_varParamStr.py --table stars_mlt_part_0870 \
--out_dir ${output_dir} \
--chunk_size 100000 --seed 88 >& stdout_0870_assign.txt &

sleep 5

nohup nice python assign_varParamStr.py --table stars_mlt_part_1100 \
--out_dir ${output_dir} \
--chunk_size 100000 --seed 112 >& stdout_1100_assign.txt &

sleep 5

nohup nice python assign_varParamStr.py --table stars_mlt_part_1160 \
--out_dir ${output_dir} \
--chunk_size 100000 --seed 321 >& stdout_1160_assign.txt &

sleep 5

nohup nice python assign_varParamStr.py --table stars_mlt_part_1180 \
--out_dir ${output_dir} \
--chunk_size 100000 --seed 425 >& stdout_1180_assign.txt &

sleep 5

nohup nice python assign_varParamStr.py --table stars_mlt_part_1200 \
--out_dir ${output_dir} \
--chunk_size 100000 --seed 5542 >& stdout_1200_assign.txt &

sleep 5

nohup nice python assign_varParamStr.py --table stars_mlt_part_1220 \
--out_dir ${output_dir} \
--chunk_size 100000 --seed 6782 >& stdout_1220_assign.txt &

sleep 5

nohup nice python assign_varParamStr.py --table stars_mlt_part_1250 \
--out_dir ${output_dir} \
--chunk_size 100000 --seed 371 >& stdout_1250_assign.txt &

sleep 5

nohup nice python assign_varParamStr.py --table stars_mlt_part_1400 \
--out_dir ${output_dir} \
--chunk_size 100000 --seed 775 >& stdout_1400_assign.txt

