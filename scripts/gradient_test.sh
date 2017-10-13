#!/bin/bash
echo "" > result.out
# vel_real.rsf -> v 
# smvel.rsf -> v_0
# dvel.rsf -> dv
# shots_d_0.rsf -> d_0
# shots_bg.rsf -> d
# alpha
# vel.rsf -> v_true
# shots_d_t -> d_t
# shots_d_b -> d_b
# err -> error
# results.out -> (x,y), using matlab to plot

# first you should generate vel_real.rsf
# sfcp < vel.rsf > vel_real.rsf
< vel_real.rsf /home/cbw/soft/madagascar-1.7/bin/sfmath output="1/input" | /home/cbw/soft/madagascar-1.7/bin/sfsmooth rect1=25 rect2=60 | /home/cbw/soft/madagascar-1.7/bin/sfmath output="1/input" > smvel.rsf 
sfmath x=vel.rsf y=smvel.rsf output="x-y" > dvel.rsf
rm -f ad_shots.rsf shots_bg.rsf
./run.sh born 8
sfcp < shots_bg.rsf > shots_d_0.rsf
sfcp < ad_shots.rsf > shots_d.rsf
for((i=-10;$i<12;i=$i+2));do
	alpha=`echo "scale=1;$i/10.0" | bc`
	echo $alpha
	out=x+$alpha*y
	sfmath x=smvel.rsf y=dvel.rsf output=$out > vel.rsf
	rm -f shots.rsf
	./run.sh fm 16
	sfmath x=shots.rsf y=shots_d_0.rsf output="x-y" > shots_d_t.rsf
	out=$alpha*x
	sfmath x=shots_d.rsf output=$out > shots_d_b.rsf
	err=`sfmath x=shots_d_t.rsf y=shots_d_b.rsf output="(x-y)*(x-y)" | sfattr want=mean | awk '{print $3}'`
	echo $alpha $err >> result.out
done
