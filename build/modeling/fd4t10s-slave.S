	.set noreorder
	.set noat
	.arch sw3
	#  /usr/sw-mpp/swcc/lib/gcc-lib/sw_64-swcc-linux/5.421-sw-485/be::5.421-sw-485

	#-----------------------------------------------------------
	# Compiling src/modeling/fd4t10s-slave.c (/tmp/ccI#.zIDAsc)
	#-----------------------------------------------------------

	#-----------------------------------------------------------
	# Options:
	#-----------------------------------------------------------
	#  Target:SW3, ISA:ISA_1, Endian:little, Pointer Size:64
	#  -O2	(Optimization level)
	#  -g0	(Debug level)
	#  -m2	(Report advisories)
	#-----------------------------------------------------------

	.file	1	"src/modeling/fd4t10s-slave.c"
	.file	2	"/usr/sw-mpp/swcc/lib/gcc-lib/sw_64-swcc-linux/5.421-sw-485/include/simd.h"
	.file	3	"/usr/sw-mpp/swcc/lib/gcc-lib/sw_64-swcc-linux/5.421-sw-485/include/dma.h"


	.section .tdata, "waT", "progbits"
	.align	0

	.section .tdata_local_fix, "waT", "progbits"
	.align	5

	.section .text1, "ax", "progbits"
	.align	4

	.section .rodata, "a", "progbits"
	.align	3
	.section .text1

	# Program Unit: fd4t10s_sw
	.align 4
	.ent	fd4t10s_sw#
	.globl	fd4t10s_sw
fd4t10s_sw:	# 0x0
	# anon36 = 32
	# anon37 = 36
	# anon38 = 40
	# anon39 = 44
	# anon40 = 224
	# anon41 = 48
	# anon42 = 52
	# anon43 = 64
	# anon44 = 96
	# anon45 = 128
	# anon46 = 160
	# anon47 = 192
	# anon48 = 256
	# anon49 = 288
	# anon54 = 320
	# return_address = 0
	# _temp_gra_spill0 = 328
	# _temp_gra_spill1 = 336
	# _temp_gra_spill2 = 344
	# _temp_gra_spill3 = 352
	# _temp_gra_spill4 = 360
	# _temp_gra_spill5 = 368
	# _temp_gra_spill6 = 376
	# _temp_gra_spill7 = 384
	# _temp_gra_spill8 = 392
	# _temp_gra_spill9 = 400
	# _temp_gra_spill10 = 408
	# _temp_gra_spill11 = 416
	# _temp_gra_spill12 = 424
	# _temp_gra_spill13 = 432
	# _temp_gra_spill14 = 440
	# _temp_gra_spill15 = 448
	# _temp_gra_spill16 = 456
	# _temp_gra_spill17 = 464
	# _temp_gra_spill18 = 472
	# _temp_gra_spill19 = 480
	# _temp_gra_spill20 = 488
	# _temp_gra_spill21 = 496
	# _temp_gra_spill22 = 504
	# _temp_gra_spill23 = 512
	# _temp_gra_spill24 = 520
	# _temp_gra_spill25 = 528
	# _temp_gra_spill26 = 536
	# _temp_gra_spill27 = 544
	# _temp_gra_spill28 = 552
	# _temp_gra_spill29 = 560
	# _temp_gra_spill30 = 568
	# _temp_gra_spill31 = 576
	# _temp_gra_spill32 = 584
	# _temp_gra_spill33 = 592
	# _temp_gra_spill34 = 600
	# _temp_gra_spill35 = 608
	# _temp_gra_spill36 = 616
	# _temp_gra_spill37 = 624
	# _temp_gra_spill38 = 632
	# _temp_gra_spill39 = 640
	# _temp_gra_spill40 = 648
	# _temp_gra_spill41 = 656
	# _temp_gra_spill42 = 664
	# _temp_gra_spill43 = 672
	# _temp_gra_spill44 = 680
	# _temp_gra_spill45 = 688
	# _temp_lgra_spill46 = 696
	# _temp_lgra_spill47 = 704
	.loc	1	12	0
#   8  #define COL 4
#   9  
#  10  const int d = 6;
#  11  
#  12  void fd4t10s_sw(struct Param *param) {
.BB1_fd4t10s_sw:
	.prologue
	ldih	$gp,0($27)               	!gpdisp!1	# [0] 0
	ldi	$gp,0($gp)                	!gpdisp!1	# [1] 0
	
$ng..fd4t10s_sw:
	ldi	$sp,-736($sp)             	# [3] 
.LCFI_fd4t10s_sw_adjustsp:
	stl	$16,392($sp)              	# [4] _temp_gra_spill8
	stl	$26,0($sp)                	# [4] return_address
.LCFI_fd4t10s_sw_store26:
	.loc	1	13	0
#  13  	int id = athread_get_id(-1);
	.frame $15,736,$26,0
	.mask 0x4008000,-736
	ldi	$16,-1($31)               	# [4] 
	.loc	1	12	0
	stl	$15,360($sp)              	# [5] _temp_gra_spill4
	stl	$14,424($sp)              	# [5] _temp_gra_spill12
	stl	$13,448($sp)              	# [6] _temp_gra_spill15
#	.spillsp $15, 0
	stl	$12,472($sp)              	# [6] _temp_gra_spill18
#	.spillsp $14, 0
	stl	$11,496($sp)              	# [7] _temp_gra_spill21
#	.spillsp $13, 0
	stl	$10,344($sp)              	# [7] _temp_gra_spill2
#	.spillsp $12, 0
	stl	$9,376($sp)               	# [8] _temp_gra_spill6
	.loc	1	13	0
#	.spillsp $11, 0
	.globl	athread_get_id
	bsr	$26,athread_get_id        	# [8] athread_get_id
.BB2_fd4t10s_sw:
	.loc	1	16	0
#  14  	int rid = id / ROW;
#  15  	int cid = id % ROW;
#  16  	Type *prev_wave = param->prev_wave;
#	.spillsp $10, 0
	ldl	$2,392($sp)               	# [0] _temp_gra_spill8
	.loc	1	38	0
#  34  	snxend = snxend < nx ? snxend : nx;
#  35  	snzend = snzend < nz ? snzend : nz;
#  36  	snx = snxend - snxbeg;
#  37  	snz = snzend - snzbeg;
#  38  	Type *local_vel = (Type *)ldm_malloc(sizeof(Type) * (snz - 2 * d) * (snx - 2 * d));
#	.spillsp $9, 0
	subw	$31,$0,$28               	# [0] 
	addw	$0,15,$6                 	# [0] 
	.loc	1	13	0
	ldih	$gp,0($26)               	!gpdisp!2	# [0] 0
	.loc	1	38	0
	selge	$0,$0,$28,$28           	# [1] 
	selge	$0,$0,$6,$6             	# [1] 
	.loc	1	13	0
	ldi	$gp,0+12($gp)             	!gpdisp!2	# [1] 0
	.loc	1	38	0
	and	$28,15,$28                	# [2] 
	sra	$6,4,$6                   	# [2] 
	.loc	1	16	0
	ldl	$1,0($2)                  	# [3] id:690
	.loc	1	20	0
	ldw	$4,32($2)                 	# [3] id:693
	.loc	1	38	0
	subw	$31,$28,$3               	# [3] 
	selge	$0,$28,$3,$28           	# [4] 
	addw	$28,1,$3                 	# [5] 
	.loc	1	16	0
	stl	$1,432($sp)               	# [6] _temp_gra_spill13
	.loc	1	17	0
	ldl	$1,8($2)                  	# [6] id:691
	.loc	1	38	0
	subw	$4,12,$5                 	# [6] 
	addw	$5,15,$8                 	# [7] 
	selge	$5,$5,$8,$5             	# [8] 
	.loc	1	17	0
	stl	$1,440($sp)               	# [9] _temp_gra_spill14
	.loc	1	18	0
	ldl	$1,16($2)                 	# [9] id:692
	.loc	1	38	0
	sra	$5,4,$5                   	# [9] 
	.loc	1	21	0
	ldw	$2,36($2)                 	# [10] id:694
	.loc	1	38	0
	addw	$5,1,$5                  	# [10] 
	mulw	$5,$3,$3                 	# [11] 
	.loc	1	18	0
	stl	$1,456($sp)               	# [12] _temp_gra_spill16
	.loc	1	38	0
	addw	$6,1,$1                  	# [12] 
	mulw	$28,$5,$28               	# [12] 
	subw	$2,12,$0                 	# [13] 
	.loc	1	21	0
	stl	$2,464($sp)               	# [13] _temp_gra_spill17
	.loc	1	38	0
	addw	$0,3,$7                  	# [14] 
	selge	$0,$0,$7,$0             	# [15] 
	sra	$0,2,$0                   	# [16] 
	addw	$0,1,$0                  	# [17] 
	mulw	$0,$1,$1                 	# [18] 
	addw	$3,12,$3                 	# [18] 
	mulw	$6,$0,$0                 	# [19] 
	cmplt	$4,$3,$16               	# [19] 
	selge	$28,$28,$31,$28         	# [19] 
	seleq	$16,$3,$4,$16           	# [20] 
	stl	$28,480($sp)              	# [20] _temp_gra_spill19
	subw	$16,$28,$16              	# [21] 
	stl	$16,584($sp)              	# [22] _temp_gra_spill32
	subw	$16,12,$16               	# [22] 
	stl	$16,616($sp)              	# [23] _temp_gra_spill36
	addw	$1,12,$1                 	# [25] 
	cmplt	$2,$1,$27               	# [26] 
	selge	$0,$0,$31,$0            	# [26] 
	seleq	$27,$1,$2,$27           	# [27] 
	mov	$0,$1                     	# [27] 
	subw	$27,$0,$27               	# [28] 
	stl	$1,488($sp)               	# [28] _temp_gra_spill20
	stl	$27,672($sp)              	# [29] _temp_gra_spill43
	subw	$27,12,$27               	# [29] 
	mull	$16,$27,$16              	# [30] 
	stl	$27,352($sp)              	# [30] _temp_gra_spill3
	stl	$27,656($sp)              	# [30] _temp_gra_spill41
	s4addl	$16,$31,$16            	# [37] 
	stl	$16,504($sp)              	# [38] _temp_gra_spill22
	.globl	ldm_malloc
	bsr	$26,ldm_malloc            	# [38] ldm_malloc
.BB3_fd4t10s_sw:
	.loc	1	39	0
#  39  	Type *local_prev_wave = (Type *)ldm_malloc(sizeof(Type) * (snz - 2 * d) * (snx - 2 * d));
	ldl	$16,504($sp)              	# [0] _temp_gra_spill22
	.loc	1	38	0
	ldih	$gp,0($26)               	!gpdisp!3	# [0] 0
	stl	$0,416($sp)               	# [1] _temp_gra_spill11
	ldi	$gp,0+4($gp)              	!gpdisp!3	# [1] 0
	.loc	1	39	0
	bsr	$26,ldm_malloc            	# [1] ldm_malloc
.BB4_fd4t10s_sw:
	.loc	1	40	0
#  40  	Type *local_curr_wave = (Type *)ldm_malloc(sizeof(Type) * snz * snx);
	ldl	$27,672($sp)              	# [0] _temp_gra_spill43
	ldl	$16,584($sp)              	# [0] _temp_gra_spill32
	.loc	1	39	0
	ldih	$gp,0($26)               	!gpdisp!4	# [0] 0
	stl	$0,400($sp)               	# [1] _temp_gra_spill9
	ldi	$gp,0+8($gp)              	!gpdisp!4	# [1] 0
	.loc	1	40	0
	mull	$16,$27,$16              	# [3] 
	stl	$27,368($sp)              	# [3] _temp_gra_spill5
	s4addl	$16,$31,$16            	# [10] 
	stl	$16,512($sp)              	# [11] _temp_gra_spill23
	bsr	$26,ldm_malloc            	# [11] ldm_malloc
.BB5_fd4t10s_sw:
	.loc	1	41	0
#  41  	Type *local_u2 = (Type *)ldm_malloc(sizeof(Type) * snz * snx);
	ldl	$16,512($sp)              	# [0] _temp_gra_spill23
	.loc	1	40	0
	ldih	$gp,0($26)               	!gpdisp!5	# [0] 0
	stl	$0,384($sp)               	# [1] _temp_gra_spill7
	ldi	$gp,0+4($gp)              	!gpdisp!5	# [1] 0
	.loc	1	41	0
	bsr	$26,ldm_malloc            	# [1] ldm_malloc
.BB6_fd4t10s_sw:
	.loc	1	73	0
#  69  
#  70  static __inline void 
#  71  dma_set_reply(dma_desc *a, unsigned int volatile *reply)
#  72  {
#  73  	((dma_desc_u *)a)->packed.reply_addr = (unsigned short)(long)reply;
	ldi	$1,48($sp)                	# [0] anon41
	.loc	1	67	0
	ldl	$16,672($sp)              	# [0] _temp_gra_spill43
	.loc	1	73	0
	ldi	$2,52($sp)                	# [0] anon42
	.loc	1	67	0
	ldi	$22,0($31)                	# [0] 
	.loc	1	73	0
	sth	$1,264($sp)               	# [1] anon48+8
	.loc	1	103	0
#  99    dma_set_bsize(&dma_put_prev_one, (snz - 2 * d) * sizeof(Type));
# 100    dma_set_stepsize(&dma_put_prev_one, (nz - (snz - 2 * d)) *sizeof(Type));
# 101  
# 102  	Type *g_curr_wave = curr_wave + snxbeg * nz + snzbeg;
# 103  	Type *g_prev_wave = prev_wave + (snxbeg + d) * nz + (snzbeg + d);
	ldl	$8,464($sp)               	# [1] _temp_gra_spill17
	.loc	1	73	0
	ldi	$1,224($sp)               	# [1] anon40
	.loc	1	67	0
	ldi	$24,-1($31)               	# [1] 
	.loc	1	73	0
	sth	$1,200($sp)               	# [2] anon47+8
	ldi	$1,44($sp)                	# [2] anon39
	.loc	1	67	0
	ldih	$22,-256($22)            	# [2] 
	ldih	$24,256($24)             	# [2] 
	.loc	1	73	0
	sth	$1,168($sp)               	# [3] anon46+8
	ldi	$1,36($sp)                	# [3] anon37
	.loc	1	79	0
#  75  
#  76  static __inline void 
#  77  dma_set_op(dma_desc *a, unsigned int op)
#  78  {
#  79  	((dma_desc_u *)a)->packed.dma_op = op;
	ldi	$7,-16($31)               	# [3] 
	.loc	1	97	0
#  93  
#  94  static __inline void 
#  95  dma_set_bsize(dma_desc *a, unsigned int bsize)
#  96  {
#  97  	((dma_desc_u *)a)->packed.stride_bsize = bsize;
	ldi	$28,-1($31)               	# [3] 
	.loc	1	73	0
	sth	$1,104($sp)               	# [4] anon44+8
	.loc	1	82	0
	ldl	$20,352($sp)              	# [4] _temp_gra_spill3
	.loc	1	76	0
	subw	$8,$16,$5                	# [4] 
	.loc	1	41	0
	mov	$0,$1                     	# [4] 
	stl	$1,408($sp)               	# [5] _temp_gra_spill10
	.loc	1	79	0
	ldbu	$3,295($sp)              	# [5] anon49+7
	.loc	1	103	0
#  99  
# 100  static __inline void 
# 101  dma_set_stepsize(dma_desc *a, unsigned int bleng)
# 102  {
# 103  	((dma_desc_u *)a)->packed.stride_length = bleng;
	addw	$5,12,$1                 	# [5] 
	.loc	1	97	0
	ldih	$28,16($28)              	# [5] 
	.loc	1	73	0
	sth	$2,296($sp)               	# [6] anon49+8
	.loc	1	67	0
	ldw	$17,128($sp)              	# [6] anon45
	.loc	1	103	0
	s4addl	$1,$31,$1              	# [6] 
	.loc	1	73	0
	ldi	$2,40($sp)                	# [6] anon38
	.loc	1	103	0
	stw	$1,300($sp)               	# [7] anon49+12
	stw	$1,268($sp)               	# [7] anon48+12
	.loc	1	82	0
	s4addl	$20,$31,$20            	# [7] 
	.loc	1	41	0
	ldih	$gp,0($26)               	!gpdisp!6	# [7] 0
	.loc	1	103	0
	stw	$1,204($sp)               	# [8] anon47+12
	stw	$1,172($sp)               	# [8] anon46+12
	s4addl	$5,$31,$5              	# [8] 
	.loc	1	41	0
	ldi	$gp,0+124($gp)            	!gpdisp!6	# [8] 0
	.loc	1	103	0
	stw	$1,140($sp)               	# [9] anon45+12
	.loc	1	73	0
	sth	$2,136($sp)               	# [9] anon45+8
	.loc	1	79	0
	ldi	$1,-16($31)               	# [9] 
	.loc	1	67	0
	and	$17,$22,$17               	# [9] 
	.loc	1	79	0
	ldbu	$2,199($sp)              	# [10] anon47+7
	.loc	1	82	0
	stl	$20,592($sp)              	# [10] _temp_gra_spill33
	.loc	1	79	0
	and	$3,$1,$3                  	# [10] 
	.loc	1	67	0
	sextw	$20,$20                 	# [10] 
	.loc	1	85	0
	zapnot	$3,1,$3                	# [11] 
	ldi	$1,-241($31)              	# [11] 
	.loc	1	67	0
	and	$20,$24,$4                	# [11] 
	.loc	1	97	0
	and	$20,$28,$20               	# [11] 
	.loc	1	85	0
	and	$3,$1,$3                  	# [12] 
	.loc	1	67	0
	or	$4,$17,$17                 	# [12] 
	.loc	1	79	0
	and	$2,$7,$2                  	# [13] 
	.loc	1	85	0
	stb	$3,295($sp)               	# [13] anon49+7
	.loc	1	97	0
	ldi	$7,0($31)                 	# [13] 
	.loc	1	67	0
	stw	$17,128($sp)              	# [14] anon45
	.loc	1	79	0
	ldbu	$17,263($sp)             	# [14] anon48+7
	.loc	1	85	0
	zapnot	$2,1,$2                	# [14] 
	.loc	1	97	0
	ldih	$7,-16($7)               	# [14] 
	.loc	1	67	0
	ldw	$18,288($sp)              	# [15] anon49
	ldw	$19,256($sp)              	# [15] anon48
	.loc	1	85	0
	and	$2,$1,$2                  	# [15] 
	.loc	1	97	0
	ldw	$3,292($sp)               	# [16] anon49+4
	.loc	1	85	0
	stb	$2,199($sp)               	# [16] anon47+7
	.loc	1	79	0
	ldi	$2,-16($31)               	# [16] 
	and	$17,$2,$17                	# [17] 
	ldbu	$2,103($sp)              	# [17] anon44+7
	.loc	1	67	0
	and	$18,$22,$18               	# [18] 
	and	$19,$22,$19               	# [18] 
	.loc	1	79	0
	or	$17,1,$17                  	# [18] 
	.loc	1	97	0
	and	$3,$7,$3                  	# [19] 
	.loc	1	67	0
	or	$4,$18,$18                 	# [19] 
	or	$4,$19,$4                  	# [19] 
	.loc	1	85	0
	zapnot	$17,1,$17              	# [19] 
	.loc	1	97	0
	or	$20,$3,$3                  	# [20] 
	.loc	1	67	0
	stw	$4,256($sp)               	# [20] anon48
	.loc	1	85	0
	and	$17,$1,$17                	# [20] 
	.loc	1	97	0
	stw	$3,292($sp)               	# [21] anon49+4
	.loc	1	79	0
	ldbu	$4,135($sp)              	# [21] anon45+7
	ldi	$3,-16($31)               	# [21] 
	and	$2,$3,$2                  	# [22] 
	ldbu	$3,167($sp)              	# [22] anon46+7
	.loc	1	85	0
	stb	$17,263($sp)              	# [23] anon48+7
	.loc	1	79	0
	ldi	$17,-16($31)              	# [23] 
	or	$2,1,$2                    	# [23] 
	and	$4,$17,$4                 	# [24] 
	ldi	$17,-16($31)              	# [24] 
	.loc	1	85	0
	zapnot	$2,1,$2                	# [24] 
	.loc	1	67	0
	stw	$18,288($sp)              	# [25] anon49
	.loc	1	97	0
	ldw	$18,196($sp)              	# [25] anon47+4
	.loc	1	79	0
	and	$3,$17,$3                 	# [25] 
	or	$4,1,$4                    	# [25] 
	or	$3,1,$3                    	# [26] 
	.loc	1	85	0
	zapnot	$4,1,$4                	# [26] 
	.loc	1	97	0
	ldw	$23,260($sp)              	# [26] anon48+4
	.loc	1	85	0
	and	$2,$1,$2                  	# [26] 
	zapnot	$3,1,$3                	# [27] 
	and	$4,$1,$4                  	# [27] 
	.loc	1	67	0
	ldl	$27,584($sp)              	# [27] _temp_gra_spill32
	.loc	1	97	0
	and	$18,$7,$18                	# [28] 
	.loc	1	85	0
	and	$3,$1,$3                  	# [28] 
	stb	$4,135($sp)               	# [29] anon45+7
	stb	$3,167($sp)               	# [29] anon46+7
	.loc	1	97	0
	or	$20,$18,$18                	# [29] 
	and	$23,$7,$23                	# [29] 
	ldw	$19,132($sp)              	# [30] anon45+4
	ldw	$17,164($sp)              	# [30] anon46+4
	.loc	1	67	0
	mulw	$27,$16,$27              	# [30] 
	.loc	1	97	0
	or	$20,$23,$23                	# [30] 
	stw	$18,196($sp)              	# [31] anon47+4
	.loc	1	67	0
	ldl	$18,368($sp)              	# [31] _temp_gra_spill5
	.loc	1	103	0
	ldl	$21,480($sp)              	# [32] _temp_gra_spill19
	.loc	1	85	0
	stb	$2,103($sp)               	# [32] anon44+7
	.loc	1	97	0
	ldw	$2,100($sp)               	# [33] anon44+4
	and	$19,$7,$19                	# [33] 
	and	$17,$7,$17                	# [33] 
	.loc	1	67	0
	ldl	$16,616($sp)              	# [34] _temp_gra_spill36
	.loc	1	97	0
	or	$20,$19,$19                	# [34] 
	or	$20,$17,$20                	# [34] 
	.loc	1	67	0
	s4addl	$18,$31,$18            	# [34] 
	ldl	$17,656($sp)              	# [35] _temp_gra_spill41
	.loc	1	102	0
	mulw	$21,$8,$26               	# [35] 
	.loc	1	67	0
	sextw	$18,$4                  	# [35] 
	.loc	1	103	0
	addw	$21,6,$21                	# [35] 
	.loc	1	79	0
	ldbu	$3,71($sp)               	# [36] anon43+7
	.loc	1	97	0
	and	$4,$28,$28                	# [36] 
	and	$2,$7,$2                  	# [36] 
	.loc	1	103	0
	mulw	$8,$21,$21               	# [36] 
	.loc	1	97	0
	or	$28,$2,$2                  	# [37] 
	stw	$19,132($sp)              	# [37] anon45+4
	.loc	1	67	0
	s4addl	$27,$31,$27            	# [37] 
	and	$4,$24,$4                 	# [37] 
	mulw	$17,$16,$17              	# [38] 
	.loc	1	97	0
	stw	$2,100($sp)               	# [38] anon44+4
	.loc	1	79	0
	ldi	$2,-16($31)               	# [38] 
	.loc	1	67	0
	sextw	$27,$27                 	# [38] 
	.loc	1	79	0
	and	$3,$2,$3                  	# [39] 
	.loc	1	102	0
	ldl	$2,488($sp)               	# [39] _temp_gra_spill20
	.loc	1	67	0
	and	$27,$24,$27               	# [39] 
	.loc	1	102	0
	ldl	$19,440($sp)              	# [40] _temp_gra_spill14
	.loc	1	79	0
	or	$3,1,$3                    	# [40] 
	.loc	1	97	0
	stw	$23,260($sp)              	# [41] anon48+4
	.loc	1	105	0
# 104  	Type *g_prev_wave_read = g_prev_wave;
# 105  	Type *g_vel = vel + (snxbeg + d) * nz + (snzbeg + d);
	ldl	$16,456($sp)              	# [41] _temp_gra_spill16
	.loc	1	85	0
	zapnot	$3,1,$3                	# [41] 
	.loc	1	103	0
	ldl	$23,432($sp)              	# [42] _temp_gra_spill13
	.loc	1	102	0
	addl	$26,$2,$26               	# [42] 
	.loc	1	85	0
	and	$3,$1,$3                  	# [42] 
	.loc	1	67	0
	ldw	$1,64($sp)                	# [43] anon43
	.loc	1	43	0
	stw	$31,32($sp)               	# [43] anon36
	.loc	1	102	0
	s4addl	$26,$19,$26            	# [43] 
	.loc	1	97	0
	stw	$20,164($sp)              	# [44] anon46+4
	.loc	1	85	0
	stb	$3,71($sp)                	# [44] anon43+7
	.loc	1	105	0
	s4addl	$21,$16,$16            	# [44] 
	.loc	1	67	0
	ldw	$25,96($sp)               	# [45] anon44
	ldw	$19,192($sp)              	# [45] anon47
	s4addl	$17,$31,$17            	# [45] 
	.loc	1	103	0
	s4addl	$21,$23,$21            	# [45] 
	.loc	1	67	0
	ldw	$20,160($sp)              	# [46] anon46
	.loc	1	97	0
	ldw	$3,68($sp)                	# [46] anon43+4
	.loc	1	105	0
	s4addl	$2,$16,$16             	# [46] 
	.loc	1	67	0
	and	$1,$22,$1                 	# [46] 
	.loc	1	43	0
	stw	$31,36($sp)               	# [47] anon37
	.loc	1	105	0
	ldl	$6,384($sp)               	# [47] _temp_gra_spill7
	.loc	1	67	0
	sextw	$17,$17                 	# [47] 
	.loc	1	103	0
	s4addl	$2,$21,$2              	# [47] 
	.loc	1	43	0
	stw	$31,40($sp)               	# [48] anon38
	.loc	1	67	0
	and	$17,$24,$17               	# [48] 
	and	$25,$22,$25               	# [48] 
	or	$1,$4,$1                   	# [48] 
	.loc	1	43	0
	stw	$31,44($sp)               	# [49] anon39
	.loc	1	67	0
	and	$20,$22,$20               	# [49] 
	.loc	1	97	0
	and	$3,$7,$3                  	# [49] 
	.loc	1	103	0
	ldi	$2,24($2)                 	# [49] 
	.loc	1	67	0
	and	$19,$22,$19               	# [50] 
	.loc	1	103	0
	stw	$5,76($sp)                	# [50] anon43+12
	.loc	1	67	0
	or	$17,$20,$20                	# [50] 
	.loc	1	105	0
	ldi	$16,24($16)               	# [50] 
	.loc	1	103	0
	stw	$5,108($sp)               	# [51] anon44+12
	.loc	1	67	0
	stw	$1,64($sp)                	# [51] anon43
	or	$25,$27,$25                	# [51] 
	.loc	1	97	0
	or	$28,$3,$28                 	# [51] 
	.loc	1	43	0
	stw	$31,48($sp)               	# [52] anon41
	.loc	1	67	0
	or	$17,$19,$17                	# [52] 
	.loc	1	73	0
	ldi	$5,32($sp)                	# [52] anon36
	.loc	1	104	0
	mov	$2,$1                     	# [52] 
	.loc	1	102	0
	stl	$26,632($sp)              	# [53] _temp_gra_spill38
	stl	$26,520($sp)              	# [53] _temp_gra_spill24
	.loc	1	105	0
	mov	$31,$7                    	# [53] 
	stl	$6,648($sp)               	# [54] _temp_gra_spill40
	.loc	1	103	0
	stl	$2,664($sp)               	# [54] _temp_gra_spill42
	.loc	1	67	0
	stw	$20,160($sp)              	# [55] anon46
	.loc	1	105	0
	stl	$16,640($sp)              	# [55] _temp_gra_spill39
	.loc	1	67	0
	stw	$25,96($sp)               	# [56] anon44
	.loc	1	97	0
	stw	$28,68($sp)               	# [56] anon43+4
	.loc	1	43	0
	stw	$31,52($sp)               	# [57] anon42
	.loc	1	67	0
	stw	$17,192($sp)              	# [57] anon47
	.loc	1	73	0
	sth	$5,72($sp)                	# [58] anon43+8
	.loc	1	104	0
	stl	$1,600($sp)               	# [58] _temp_gra_spill34
	.align	4
.Lt_0_16130:
#<loop> Loop body line 105, nesting depth: 1, iterations: 11
	.loc	1	112	0
# 108  	//dma(dma_get_prev_wave, (long)g_prev_wave, (long)local_prev_wave);
# 109  	//dma_wait(&replyget_prev_wave, 1);
# 110  
# 111  	for(ix = 0 ; ix < 2 * d - 2 + 1 ; ix ++) {
# 112  		dma(dma_get_curr_one, (long)(g_curr_wave), (long)(local_curr_wave + ix * snz));
	ldw	$1,32($sp)                	# [0] anon36
	vldd	$2,64($sp)               	# [0] anon43
	ldl	$3,632($sp)               	# [1] _temp_gra_spill38
	cmpeq	$1,1,$1                 	# [3] 
	dma      $2,$3,$6             	# [4] 
	bne	$1,.Lt_0_16386            	# [4] 
.BB49_fd4t10s_sw:
.BB47_fd4t10s_sw:
#<loop> Loop body line 112
#<loop> unrolled 3 times
	ldw	$1,32($sp)                	# [0] anon36
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.BB48_fd4t10s_sw       	# [4] 
.BB50_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .BB49_fd4t10s_sw
#<loop> unrolled 3 times
	ldw	$1,32($sp)                	# [0] anon36
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.BB48_fd4t10s_sw       	# [4] 
.Lt_0_16898:
#<loop> Part of loop body line 65536, head labeled .BB49_fd4t10s_sw
	ldw	$1,32($sp)                	# [0] anon36
	cmpeq	$1,1,$1                 	# [3] 
	beq	$1,.BB49_fd4t10s_sw       	# [4] 
.Lt_0_16386:
.BB48_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_16130
	.loc	1	115	0
# 113  		dma_wait(&replyget_curr_one, 1);
# 114  		replyget_curr_one = 0;
# 115  		g_curr_wave += nz;
	ldl	$16,632($sp)              	# [0] _temp_gra_spill38
	.loc	1	111	0
	addw	$7,1,$7                  	# [0] 
	.loc	1	115	0
	s4addl	$8,$31,$2              	# [0] 
	.loc	1	111	0
	addl	$6,$18,$6                	# [0] 
	cmpeq	$7,11,$1                	# [1] 
	.loc	1	114	0
	stw	$31,32($sp)               	# [1] anon36
	.loc	1	115	0
	stl	$2,624($sp)               	# [2] _temp_gra_spill37
	s4addl	$8,$16,$16             	# [3] 
	stl	$16,632($sp)              	# [4] _temp_gra_spill38
	.loc	1	111	0
	beq	$1,.Lt_0_16130            	# [4] 
.BB10_fd4t10s_sw:
	.loc	1	132	0
# 128  	dma_wait(&replyget_curr_wave, 1);
# 129  	replyget_curr_wave = 0;
# 130  	*/
# 131  
# 132    Type *a = (Type *)ldm_malloc(sizeof(Type) * 6);
	mov	24,$16                    	# [0] 
	bsr	$26,ldm_malloc            	# [0] ldm_malloc
.BB11_fd4t10s_sw:
	.loc	1	138	0
# 134    a[1] = +1.78858721;
# 135    a[2] = -0.31660756;
# 136    a[3] = +0.07612173;
# 137    a[4] = -0.01626042;
# 138    a[5] = +0.00216736;
	ldl	$5,584($sp)               	# [0] _temp_gra_spill32
	.loc	1	133	0
	ldi	$24,-29876($31)           	# [0] 
	.loc	1	134	0
	ldi	$22,-25074($31)           	# [0] 
	.loc	1	135	0
	ldi	$19,17228($31)            	# [0] 
	.loc	1	136	0
	ldi	$21,31927($31)            	# [1] 
	.loc	1	137	0
	ldi	$23,-22894($31)           	# [1] 
	.loc	1	138	0
	ldi	$20,-16055($31)           	# [1] 
	.loc	1	133	0
	ldi	$17,0($31)                	# [1] 
	.loc	1	138	0
	ldl	$25,408($sp)              	# [2] _temp_gra_spill10
	.loc	1	134	0
	ldi	$8,0($31)                 	# [2] 
	.loc	1	135	0
	ldi	$16,0($31)                	# [2] 
	.loc	1	136	0
	ldi	$7,0($31)                 	# [2] 
	.loc	1	138	0
	ldi	$2,0($31)                 	# [3] 
	.loc	1	133	0
	ldih	$24,16377($24)           	# [3] 
	.loc	1	134	0
	ldih	$22,16381($22)           	# [3] 
	.loc	1	135	0
	ldih	$19,-16428($19)          	# [3] 
	.loc	1	136	0
	ldih	$21,16307($21)           	# [4] 
	.loc	1	137	0
	ldih	$23,-16495($23)          	# [4] 
	.loc	1	138	0
	ldih	$20,16226($20)           	# [4] 
	.loc	1	135	0
	sll	$19,32,$19                	# [4] 
	.loc	1	133	0
	sll	$24,32,$24                	# [5] 
	.loc	1	134	0
	sll	$22,32,$22                	# [5] 
	.loc	1	135	0
	ldih	$16,24576($16)           	# [5] 
	.loc	1	138	0
	subw	$5,5,$6                  	# [5] 
	.loc	1	136	0
	sll	$21,32,$21                	# [6] 
	.loc	1	137	0
	sll	$23,32,$23                	# [6] 
	.loc	1	133	0
	ldih	$17,-8192($17)           	# [6] 
	.loc	1	134	0
	ldih	$8,-24576($8)            	# [6] 
	.loc	1	138	0
	sll	$20,32,$20                	# [7] 
	.loc	1	136	0
	ldih	$7,-16384($7)            	# [7] 
	.loc	1	138	0
	ldih	$2,-32768($2)            	# [7] 
	.loc	1	132	0
	ldih	$gp,0($26)               	!gpdisp!7	# [7] 0
	.loc	1	233	0
# 229                             local_curr_wave[curPos - 5 * snz]  +  local_curr_wave[curPos + 5 * snz]);
# 230  		}
# 231  	}
# 232  
# 233  	dma(dma_put_prev_one, (long)g_prev_wave_read, (long)(local_prev_wave)); //just for start
	subw	$5,6,$5                  	# [8] 
	.loc	1	135	0
	addl	$16,$19,$19              	# [8] 
	.loc	1	132	0
	mov	$0,$18                    	# [8] 
	.loc	1	138	0
	cmple	$6,5,$1                 	# [8] 
	.loc	1	133	0
	addl	$17,$24,$17              	# [9] 
	.loc	1	134	0
	addl	$8,$22,$8                	# [9] 
	.loc	1	138	0
	addl	$2,$20,$2                	# [9] 
	.loc	1	136	0
	addl	$7,$21,$7                	# [9] 
	.loc	1	137	0
	addl	$16,$23,$16              	# [10] 
	.loc	1	138	0
	stl	$25,608($sp)              	# [10] _temp_gra_spill35
	.loc	1	132	0
	ldi	$gp,0+124($gp)            	!gpdisp!7	# [10] 0
	.loc	1	233	0
	cmple	$5,6,$5                 	# [10] 
	.loc	1	135	0
	fsts	$19,8($0)                	# [11] id:775
	.loc	1	133	0
	stl	$18,680($sp)              	# [11] _temp_gra_spill44
	fsts	$17,0($0)                	# [12] id:773
	.loc	1	134	0
	fsts	$8,4($0)                 	# [12] id:774
	.loc	1	138	0
	fsts	$2,20($0)                	# [13] id:778
	.loc	1	136	0
	fsts	$7,12($0)                	# [13] id:776
	.loc	1	137	0
	fsts	$16,16($0)               	# [14] id:777
	.loc	1	138	0
	bne	$1,.Lt_0_26882            	# [14] 
.BB12_fd4t10s_sw:
	ldl	$0,416($sp)               	# [0] _temp_gra_spill11
	ldl	$10,672($sp)              	# [0] _temp_gra_spill43
	ldi	$28,3073($31)             	# [0] 
	mov	5,$17                     	# [0] 
	ldl	$16,400($sp)              	# [1] _temp_gra_spill9
	addw	$6,6,$24                 	# [1] 
	sll	$28,52,$28                	# [1] 
	mov	$31,$15                   	# [1] 
	stl	$5,328($sp)               	# [2] _temp_gra_spill0
	stl	$24,576($sp)              	# [2] _temp_gra_spill31
	mov	11,$14                    	# [2] 
	stl	$0,688($sp)               	# [3] _temp_gra_spill45
	ldl	$0,584($sp)               	# [3] _temp_gra_spill32
	s4addl	$10,$10,$18            	# [3] 
	addl	$10,$10,$20              	# [3] 
	s4subl	$10,$10,$21            	# [4] 
	sextw	$18,$18                 	# [4] 
	subw	$10,5,$19                	# [4] 
	s4addl	$10,$31,$23            	# [4] 
	stl	$18,536($sp)              	# [5] _temp_gra_spill26
	addw	$18,5,$13                	# [5] 
	subw	$10,10,$25               	# [5] 
	sextw	$20,$20                 	# [5] 
	subw	$0,10,$22                	# [6] 
	sextw	$21,$21                 	# [6] 
	sextw	$23,$23                 	# [6] 
	addw	$19,$18,$18              	# [6] 
	stl	$25,544($sp)              	# [7] _temp_gra_spill27
	stl	$20,568($sp)              	# [7] _temp_gra_spill30
	addw	$10,$13,$11              	# [7] 
	subw	$13,$10,$12              	# [7] 
	stl	$22,336($sp)              	# [8] _temp_gra_spill1
	stl	$21,560($sp)              	# [8] _temp_gra_spill29
	cmplt	$17,$19,$17             	# [8] 
	subw	$0,13,$0                 	# [8] 
	stl	$23,552($sp)              	# [9] _temp_gra_spill28
	addw	$10,$18,$10              	# [9] 
	br	$31,.Lt_0_18178            	# [9] 
	.align	4
.Lt_0_18434:
#<loop> Part of loop body line 65536, head labeled .Lt_0_18178
	.loc	1	151	0
	ldl	$5,520($sp)               	# [0] _temp_gra_spill24
.Lt_0_26370:
#<loop> Part of loop body line 65536, head labeled .Lt_0_18178
	.loc	1	152	0
	ldw	$1,32($sp)                	# [0] anon36
	.loc	1	151	0
	mov	$5,$2                     	# [0] 
	stl	$2,632($sp)               	# [1] _temp_gra_spill38
	.loc	1	152	0
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.Lt_0_18690            	# [4] 
.BB53_fd4t10s_sw:
.BB51_fd4t10s_sw:
#<loop> Loop body line 152
#<loop> unrolled 3 times
	ldw	$1,32($sp)                	# [0] anon36
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.BB52_fd4t10s_sw       	# [4] 
.BB54_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .BB53_fd4t10s_sw
#<loop> unrolled 3 times
	ldw	$1,32($sp)                	# [0] anon36
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.BB52_fd4t10s_sw       	# [4] 
.Lt_0_19202:
#<loop> Part of loop body line 65536, head labeled .BB53_fd4t10s_sw
	ldw	$1,32($sp)                	# [0] anon36
	cmpeq	$1,1,$1                 	# [3] 
	beq	$1,.BB53_fd4t10s_sw       	# [4] 
.Lt_0_18690:
.BB52_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_18178
	.loc	1	153	0
	stw	$31,32($sp)               	# [0] anon36
	.loc	1	159	0
	ldw	$1,40($sp)                	# [0] anon38
	.loc	1	162	0
	ldl	$5,688($sp)               	# [1] _temp_gra_spill45
	vldd	$2,256($sp)              	# [1] anon48
	ldl	$3,640($sp)               	# [2] _temp_gra_spill39
	.loc	1	159	0
	cmpeq	$1,1,$1                 	# [3] 
	.loc	1	162	0
	s4addl	$7,$5,$4               	# [4] 
	dma      $2,$3,$4             	# [5] 
	.loc	1	159	0
	bne	$1,.Lt_0_19714            	# [5] 
.BB57_fd4t10s_sw:
.BB55_fd4t10s_sw:
#<loop> Loop body line 159
#<loop> unrolled 3 times
	ldw	$1,40($sp)                	# [0] anon38
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.BB56_fd4t10s_sw       	# [4] 
.BB58_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .BB57_fd4t10s_sw
#<loop> unrolled 3 times
	ldw	$1,40($sp)                	# [0] anon38
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.BB56_fd4t10s_sw       	# [4] 
.Lt_0_20226:
#<loop> Part of loop body line 65536, head labeled .BB57_fd4t10s_sw
	ldw	$1,40($sp)                	# [0] anon38
	cmpeq	$1,1,$1                 	# [3] 
	beq	$1,.BB57_fd4t10s_sw       	# [4] 
.Lt_0_19714:
.BB56_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_18178
	.loc	1	160	0
	stw	$31,40($sp)               	# [0] anon38
	.loc	1	164	0
	ldw	$1,48($sp)                	# [0] anon41
	.loc	1	162	0
	mov	$6,$2                     	# [0] 
	stl	$2,640($sp)               	# [1] _temp_gra_spill39
	.loc	1	164	0
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.Lt_0_20738            	# [4] 
.BB61_fd4t10s_sw:
.BB59_fd4t10s_sw:
#<loop> Loop body line 164
#<loop> unrolled 3 times
	ldw	$1,48($sp)                	# [0] anon41
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.BB60_fd4t10s_sw       	# [4] 
.BB62_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .BB61_fd4t10s_sw
#<loop> unrolled 3 times
	ldw	$1,48($sp)                	# [0] anon41
	cmpeq	$1,1,$1                 	# [3] 
	bne	$1,.BB60_fd4t10s_sw       	# [4] 
.Lt_0_21250:
#<loop> Part of loop body line 65536, head labeled .BB61_fd4t10s_sw
	ldw	$1,48($sp)                	# [0] anon41
	cmpeq	$1,1,$1                 	# [3] 
	beq	$1,.BB61_fd4t10s_sw       	# [4] 
.Lt_0_20738:
.BB60_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_18178
	.loc	1	165	0
	stw	$31,48($sp)               	# [0] anon41
	beq	$17,.Lt_0_21762           	# [0] 
.BB24_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_18178
	ldl	$24,536($sp)              	# [0] _temp_gra_spill26
	ldl	$21,552($sp)              	# [0] _temp_gra_spill28
	mov	$11,$7                    	# [0] 
	ldl	$20,560($sp)              	# [1] _temp_gra_spill29
	ldl	$22,568($sp)              	# [1] _temp_gra_spill30
	ldl	$1,648($sp)               	# [2] _temp_gra_spill40
	ldl	$6,608($sp)               	# [2] _temp_gra_spill35
	ldl	$2,544($sp)               	# [3] _temp_gra_spill27
	addw	$13,$24,$18              	# [3] 
	addw	$21,$13,$27              	# [3] 
	subw	$13,$24,$24              	# [3] 
	addw	$20,$13,$25              	# [4] 
	addw	$22,$13,$19              	# [4] 
	subw	$13,$21,$21              	# [4] 
	subw	$13,$20,$20              	# [4] 
	subw	$13,$22,$22              	# [5] 
	s4addl	$12,$1,$23             	# [5] 
	s4addl	$11,$1,$5              	# [5] 
	s4addl	$13,$1,$8              	# [5] 
	stl	$2,528($sp)               	# [6] _temp_gra_spill25
	s4addl	$13,$6,$6              	# [6] 
	s4addl	$18,$1,$18             	# [6] 
	s4addl	$27,$1,$27             	# [6] 
	s4addl	$24,$1,$24             	# [7] 
	s4addl	$25,$1,$25             	# [7] 
	s4addl	$19,$1,$19             	# [7] 
	s4addl	$21,$1,$21             	# [7] 
	s4addl	$20,$1,$20             	# [8] 
	s4addl	$22,$1,$22             	# [8] 
	.align	4
.Lt_0_22274:
.BB63_fd4t10s_sw:
#<loop> Loop body line 165, nesting depth: 2, estimated iterations: 100
	stl	$0,696($sp)               	# [0] _temp_lgra_spill46
	.loc	1	219	0
	flds	$2,4($8)                 	# [0] id:809
	.loc	1	167	0
	addw	$7,1,$7                  	# [0] 
	.loc	1	219	0
	flds	$0,-4($8)                	# [1] id:808
	stl	$11,704($sp)              	# [1] _temp_lgra_spill47
	fadds	$0,$2,$1                	# [5] 
	flds	$0,0($23)                	# [5] id:807
	.loc	1	167	0
	ldi	$23,4($23)                	# [5] 
	.loc	1	219	0
	fadds	$0,$1,$2                	# [9] 
	flds	$0,0($5)                 	# [9] id:806
	ldl	$1,680($sp)               	# [9] _temp_gra_spill44
	.loc	1	167	0
	ldi	$5,4($5)                  	# [9] 
	.loc	1	219	0
	flds	$3,4($1)                 	# [12] id:805
	fadds	$0,$2,$4                	# [13] 
	flds	$0,0($1)                 	# [13] id:811
	fmuls	$3,$4,$2                	# [17] 
	vfmov	$0,$4                   	# [17] 
	flds	$3,0($8)                 	# [17] id:810
	vfmov	$3,$0                   	# [21] 
	fmuld	$4,$28,$3               	# [21] 
	vfmov	$2,$4                   	# [22] 
	fmad	$0,$3,$4,$2              	# [26] 
	flds	$3,20($8)                	# [26] id:789
	flds	$0,-20($8)               	# [26] id:788
	fadds	$0,$3,$4                	# [30] 
	flds	$0,0($24)                	# [30] id:787
	.loc	1	167	0
	ldi	$24,4($24)                	# [30] 
	.loc	1	219	0
	fadds	$0,$4,$3                	# [34] 
	flds	$0,0($18)                	# [34] id:786
	.loc	1	167	0
	ldi	$18,4($18)                	# [34] 
	.loc	1	219	0
	fadds	$0,$3,$4                	# [38] 
	flds	$3,16($8)                	# [38] id:794
	flds	$0,-16($8)               	# [38] id:793
	fadds	$0,$3,$9                	# [42] 
	flds	$0,0($21)                	# [42] id:792
	.loc	1	167	0
	ldi	$21,4($21)                	# [42] 
	.loc	1	219	0
	fadds	$0,$9,$3                	# [46] 
	flds	$0,0($27)                	# [46] id:791
	.loc	1	167	0
	ldi	$27,4($27)                	# [46] 
	.loc	1	219	0
	fadds	$0,$3,$11               	# [50] 
	flds	$3,12($8)                	# [50] id:799
	flds	$0,-12($8)               	# [50] id:798
	fadds	$0,$3,$9                	# [54] 
	flds	$0,0($20)                	# [54] id:797
	.loc	1	167	0
	ldi	$20,4($20)                	# [54] 
	.loc	1	219	0
	fadds	$0,$9,$3                	# [58] 
	flds	$0,0($25)                	# [58] id:796
	flds	$9,8($8)                 	# [58] id:804
	.loc	1	167	0
	ldi	$25,4($25)                	# [58] 
	.loc	1	219	0
	fadds	$0,$3,$28               	# [62] 
	flds	$0,-8($8)                	# [62] id:803
	.loc	1	167	0
	ldi	$8,4($8)                  	# [62] 
	.loc	1	219	0
	fadds	$0,$9,$3                	# [66] 
	flds	$0,0($22)                	# [66] id:802
	.loc	1	167	0
	ldi	$22,4($22)                	# [66] 
	.loc	1	219	0
	fadds	$0,$3,$9                	# [70] 
	flds	$0,0($19)                	# [70] id:801
	.loc	1	167	0
	ldi	$19,4($19)                	# [70] 
	.loc	1	219	0
	fadds	$0,$9,$3                	# [74] 
	flds	$0,8($1)                 	# [74] id:800
	fmuls	$0,$3,$9                	# [78] 
	flds	$0,12($1)                	# [78] id:795
	fmuls	$0,$28,$3               	# [82] 
	vfmov	$9,$0                   	# [82] 
	ldi	$28,3073($31)             	# [82] 
	sll	$28,52,$28                	# [83] 
	faddd	$0,$2,$9                	# [86] 
	flds	$2,16($1)                	# [86] id:790
	fmuls	$2,$11,$0               	# [90] 
	vfmov	$3,$2                   	# [90] 
	ldl	$11,704($sp)              	# [90] _temp_lgra_spill47
	faddd	$2,$9,$3                	# [94] 
	flds	$2,20($1)                	# [94] id:785
	fmuls	$2,$4,$1                	# [98] 
	vfmov	$0,$2                   	# [98] 
	vfmov	$1,$0                   	# [102] 
	faddd	$2,$3,$1                	# [103] 
	faddd	$0,$1,$2                	# [107] 
	.loc	1	167	0
	cmpeq	$7,$10,$1               	# [107] 
	.loc	1	219	0
	fcvtds	$2,$0                  	# [111] 
	fsts	$0,0($6)                 	# [112] id:812
	ldl	$0,696($sp)               	# [112] _temp_lgra_spill46
	.loc	1	167	0
	ldi	$6,4($6)                  	# [112] 
	beq	$1,.Lt_0_22274            	# [112] 
.Lt_0_21762:
.BB64_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_18178
	.loc	1	148	0
	ldl	$18,672($sp)              	# [0] _temp_gra_spill43
	addw	$18,$11,$11              	# [3] 
	addw	$18,$12,$12              	# [3] 
	addw	$18,$13,$13              	# [3] 
	addw	$18,$10,$10              	# [3] 
	bne	$26,.Lt_0_17666           	# [4] 
.Lt_0_18178:
#<loop> Loop body line 138, nesting depth: 1, estimated iterations: 100
	.loc	1	157	0
	ldl	$8,616($sp)               	# [0] _temp_gra_spill36
	.loc	1	150	0
	ldl	$1,584($sp)               	# [0] _temp_gra_spill32
	.loc	1	157	0
	ldl	$7,656($sp)               	# [1] _temp_gra_spill41
	.loc	1	150	0
	ldl	$4,672($sp)               	# [1] _temp_gra_spill43
	.loc	1	158	0
	ldl	$5,600($sp)               	# [2] _temp_gra_spill34
	ldl	$6,624($sp)               	# [2] _temp_gra_spill37
	.loc	1	157	0
	cmplt	$15,$8,$8               	# [3] 
	.loc	1	150	0
	cmplt	$14,$1,$23              	# [3] 
	.loc	1	157	0
	seleq	$8,$0,$15,$8            	# [4] 
	.loc	1	150	0
	seleq	$23,$31,$14,$23         	# [4] 
	.loc	1	162	0
	ldl	$22,640($sp)              	# [4] _temp_gra_spill39
	.loc	1	157	0
	vldd	$18,128($sp)             	# [4] anon45
	mulw	$7,$8,$7                 	# [5] 
	.loc	1	150	0
	ldl	$20,648($sp)              	# [5] _temp_gra_spill40
	vldd	$2,64($sp)               	# [5] anon43
	.loc	1	157	0
	cmpeq	$0,$8,$8                	# [5] 
	.loc	1	150	0
	mulw	$4,$23,$4                	# [6] 
	.loc	1	148	0
	ldl	$26,576($sp)              	# [6] _temp_gra_spill31
	.loc	1	150	0
	ldl	$3,632($sp)               	# [6] _temp_gra_spill38
	.loc	1	158	0
	addl	$6,$5,$21                	# [6] 
	.loc	1	157	0
	cmpeq	$8,0,$8                 	# [7] 
	.loc	1	162	0
	addl	$6,$22,$6                	# [7] 
	.loc	1	148	0
	addw	$14,1,$14                	# [7] 
	.loc	1	151	0
	cmple	$1,$23,$1               	# [7] 
	.loc	1	162	0
	selne	$8,$6,$22,$6            	# [8] 
	.loc	1	158	0
	selne	$8,$21,$5,$8            	# [8] 
	.loc	1	148	0
	addw	$15,1,$15                	# [8] 
	cmpeq	$14,$26,$26             	# [9] 
	.loc	1	157	0
	s4addl	$7,$16,$19             	# [12] 
	dma      $18,$5,$19           	# [13] 
	.loc	1	150	0
	s4addl	$4,$20,$4              	# [13] 
	.loc	1	158	0
	mov	$8,$5                     	# [13] 
	stl	$5,600($sp)               	# [14] _temp_gra_spill34
	.loc	1	150	0
	dma      $2,$3,$4             	# [14] 
	.loc	1	151	0
	bne	$1,.Lt_0_18434            	# [14] 
.BB15_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_18178
	ldl	$5,624($sp)               	# [0] _temp_gra_spill37
	addl	$5,$3,$5                 	# [3] 
	br	$31,.Lt_0_26370            	# [3] 
.Lt_0_26882:
	.loc	1	148	0
	ldl	$0,416($sp)               	# [0] _temp_gra_spill11
	ldl	$16,400($sp)              	# [0] _temp_gra_spill9
	stl	$5,328($sp)               	# [1] _temp_gra_spill0
	stl	$0,688($sp)               	# [3] _temp_gra_spill45
.Lt_0_17666:
	.loc	1	233	0
	vldd	$15,288($sp)             	# [0] anon49
	ldl	$1,600($sp)               	# [0] _temp_gra_spill34
	ldl	$0,328($sp)               	# [1] _temp_gra_spill0
	dma      $15,$1,$16           	# [4] 
	bne	$0,.Lt_0_23042            	# [4] 
.BB31_fd4t10s_sw:
	ldl	$10,672($sp)              	# [0] _temp_gra_spill43
	ldi	$0,21845($31)             	# [0] 
	ldi	$5,0($31)                 	# [0] 
	mov	6,$13                     	# [0] 
	fcmpeq	$31,$31,$25            	# [0] 
	ldih	$0,16309($0)             	# [1] 
	ldi	$22,3073($31)             	# [1] 
	ldih	$5,24576($5)             	# [1] 
	mov	$31,$27                   	# [1] 
	sll	$0,32,$0                  	# [2] 
	sll	$22,52,$22                	# [2] 
	mov	6,$28                     	# [2] 
	mov	$16,$11                   	# [2] 
	s4subl	$10,$10,$9             	# [3] 
	subw	$10,6,$1                 	# [3] 
	addl	$5,$0,$5                 	# [3] 
	addl	$9,$9,$9                 	# [4] 
	cmplt	$13,$1,$13              	# [4] 
	sextw	$9,$9                   	# [5] 
	addw	$9,$1,$24                	# [6] 
	addw	$9,6,$9                  	# [6] 
	addw	$10,$24,$24              	# [7] 
	addw	$10,$9,$10               	# [7] 
	.align	4
.Lt_0_23554:
#<loop> Loop body line 233, nesting depth: 1, estimated iterations: 100
	.loc	1	249	0
# 245  		dma_wait(&replyput_prev_one, 1);
# 246  		replyput_prev_one = 0;
# 247  		//int ix_put_prev = ix != d ?  ix - d : 0;
# 248  		int ix_put_prev = ix - d;
# 249  		dma(dma_put_prev_one, (long)g_prev_wave, (long)(local_prev_wave + ix_put_prev * (snz - 2 * d)));
	ldl	$18,664($sp)              	# [0] _temp_gra_spill42
	.loc	1	234	0
	ldl	$17,592($sp)              	# [0] _temp_gra_spill33
	addw	$27,1,$27                	# [0] 
	.loc	1	249	0
	dma      $15,$18,$11          	# [3] 
	.loc	1	234	0
	addl	$11,$17,$11              	# [3] 
	.loc	1	233	0
	beq	$13,.BB76_fd4t10s_sw      	# [3] 
.BB34_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_23554
	ldl	$20,672($sp)              	# [0] _temp_gra_spill43
	ldl	$23,608($sp)              	# [0] _temp_gra_spill35
	subw	$28,6,$21                	# [0] 
	mov	$10,$7                    	# [0] 
	.loc	1	234	0
	ldl	$12,616($sp)              	# [1] _temp_gra_spill36
	.loc	1	233	0
	ldl	$8,648($sp)               	# [1] _temp_gra_spill40
	s4addl	$21,$16,$18            	# [1] 
	subw	$9,$20,$20               	# [3] 
	s4addl	$10,$23,$19            	# [3] 
	s4addl	$9,$23,$6              	# [3] 
	.loc	1	234	0
	cmpeq	$27,$12,$12             	# [4] 
	.loc	1	233	0
	s4addl	$9,$8,$8               	# [4] 
	s4addl	$20,$23,$20            	# [4] 
	.align	4
.Lt_0_24322:
.BB65_fd4t10s_sw:
#<loop> Loop body line 233, nesting depth: 2, estimated iterations: 100
	.loc	1	238	0
	ldl	$0,688($sp)               	# [0] _temp_gra_spill45
	.loc	1	239	0
	flds	$17,0($8)                	# [0] id:819
	.loc	1	235	0
	addw	$7,1,$7                  	# [0] 
	ldi	$8,4($8)                  	# [0] 
	.loc	1	239	0
	flds	$3,0($18)                	# [1] id:820
	.loc	1	238	0
	s4addl	$21,$0,$0              	# [3] 
	.loc	1	235	0
	addw	$21,1,$21                	# [3] 
	.loc	1	238	0
	flds	$1,0($0)                 	# [4] id:813
	.loc	1	239	0
	ldi	$0,1023($31)              	# [4] 
	vfmov	$17,$4                  	# [4] 
	sll	$0,52,$0                  	# [5] 
	vfmov	$3,$23                  	# [5] 
	fdivs	$0,$1,$2                	# [8] 
	flds	$1,4($6)                 	# [8] id:818
	flds	$0,-4($6)                	# [8] id:817
	fmsd	$4,$25,$23,$17           	# [9] 
	fadds	$0,$1,$3                	# [12] 
	flds	$1,0($6)                 	# [12] id:814
	flds	$0,0($20)                	# [12] id:816
	.loc	1	235	0
	ldi	$6,4($6)                  	# [12] 
	ldi	$20,4($20)                	# [12] 
	.loc	1	239	0
	fadds	$0,$3,$23               	# [16] 
	fmuls	$2,$1,$3                	# [16] 
	flds	$0,0($19)                	# [16] id:815
	.loc	1	235	0
	ldi	$19,4($19)                	# [16] 
	.loc	1	239	0
	fadds	$0,$23,$4               	# [20] 
	fmuls	$2,$5,$23               	# [20] 
	vfmov	$3,$0                   	# [21] 
	fmuls	$2,$23,$3               	# [24] 
	faddd	$0,$17,$2               	# [25] 
	fmas	$1,$22,$4,$0             	# [25] 
	fmuls	$0,$3,$1                	# [29] 
	vfmov	$1,$0                   	# [33] 
	faddd	$0,$2,$1                	# [37] 
	.loc	1	235	0
	cmpeq	$7,$24,$0               	# [37] 
	.loc	1	239	0
	fcvtds	$1,$23                 	# [41] 
	fsts	$23,0($18)               	# [42] id:821
	.loc	1	235	0
	ldi	$18,4($18)                	# [42] 
	beq	$0,.Lt_0_24322            	# [42] 
.Lt_0_23810:
.BB66_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_23554
	ldw	$0,52($sp)                	# [0] anon42
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.Lt_0_24834            	# [4] 
.BB69_fd4t10s_sw:
.BB67_fd4t10s_sw:
#<loop> Loop body line 235
#<loop> unrolled 3 times
	ldw	$0,52($sp)                	# [0] anon42
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.BB68_fd4t10s_sw       	# [4] 
.BB70_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .BB69_fd4t10s_sw
#<loop> unrolled 3 times
	ldw	$0,52($sp)                	# [0] anon42
	cmpeq	$0,1,$0                 	# [3] 
	bne	$0,.BB68_fd4t10s_sw       	# [4] 
.Lt_0_25346:
#<loop> Part of loop body line 65536, head labeled .BB69_fd4t10s_sw
	ldw	$0,52($sp)                	# [0] anon42
	cmpeq	$0,1,$0                 	# [3] 
	beq	$0,.BB69_fd4t10s_sw       	# [4] 
.Lt_0_24834:
.BB68_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_23554
	.loc	1	249	0
	ldl	$0,664($sp)               	# [0] _temp_gra_spill42
	ldl	$1,624($sp)               	# [0] _temp_gra_spill37
	.loc	1	246	0
	stw	$31,52($sp)               	# [1] anon42
	.loc	1	249	0
	addl	$1,$0,$0                 	# [3] 
	stl	$0,664($sp)               	# [4] _temp_gra_spill42
	.loc	1	234	0
	ldl	$0,672($sp)               	# [4] _temp_gra_spill43
	addw	$0,$10,$10               	# [7] 
	addw	$0,$9,$9                 	# [7] 
	addw	$0,$24,$24               	# [7] 
	ldl	$0,656($sp)               	# [8] _temp_gra_spill41
	addw	$0,$28,$28               	# [11] 
	beq	$12,.Lt_0_23554           	# [11] 
.Lt_0_23042:
	.loc	1	263	0
# 259  	dma_wait(&replyput_prev_wave, 1);
# 260  	replyput_prev_wave = 0;
# 261  	*/
# 262  
# 263  	ldm_free(local_prev_wave, sizeof(Type) * (snx - 2 * d) * (snz - 2 * d));
	ldl	$17,504($sp)              	# [0] _temp_gra_spill22
	.globl	ldm_free
	bsr	$26,ldm_free              	# [0] ldm_free
.BB41_fd4t10s_sw:
	.loc	1	264	0
# 264  	ldm_free(local_curr_wave, sizeof(Type) * snx * snz);
	ldl	$16,648($sp)              	# [0] _temp_gra_spill40
	ldl	$17,512($sp)              	# [0] _temp_gra_spill23
	.loc	1	263	0
	ldih	$gp,0($26)               	!gpdisp!8	# [0] 0
	ldi	$gp,0+8($gp)              	!gpdisp!8	# [1] 0
	.loc	1	264	0
	bsr	$26,ldm_free              	# [1] ldm_free
.BB42_fd4t10s_sw:
	.loc	1	265	0
# 265  	ldm_free(local_u2, sizeof(Type) * snx * snz);
	ldl	$16,608($sp)              	# [0] _temp_gra_spill35
	ldl	$17,512($sp)              	# [0] _temp_gra_spill23
	.loc	1	264	0
	ldih	$gp,0($26)               	!gpdisp!9	# [0] 0
	ldi	$gp,0+8($gp)              	!gpdisp!9	# [1] 0
	.loc	1	265	0
	bsr	$26,ldm_free              	# [1] ldm_free
.BB43_fd4t10s_sw:
	.loc	1	266	0
# 266  	ldm_free(local_vel, sizeof(Type) * (snx - 2 * d) * (snz - 2 * d));
	ldl	$16,688($sp)              	# [0] _temp_gra_spill45
	ldl	$17,504($sp)              	# [0] _temp_gra_spill22
	.loc	1	265	0
	ldih	$gp,0($26)               	!gpdisp!10	# [0] 0
	ldi	$gp,0+8($gp)              	!gpdisp!10	# [1] 0
	.loc	1	266	0
	bsr	$26,ldm_free              	# [1] ldm_free
.BB44_fd4t10s_sw:
	.loc	1	267	0
# 267  	ldm_free(a, sizeof(Type) * 6);
	ldl	$16,680($sp)              	# [0] _temp_gra_spill44
	.loc	1	266	0
	ldih	$gp,0($26)               	!gpdisp!11	# [0] 0
	.loc	1	267	0
	mov	24,$17                    	# [0] 
	.loc	1	266	0
	ldi	$gp,0+4($gp)              	!gpdisp!11	# [1] 0
	.loc	1	267	0
	bsr	$26,ldm_free              	# [1] ldm_free
.BB45_fd4t10s_sw:
	ldih	$gp,0($26)               	!gpdisp!12	# [0] 0
	ldl	$26,0($sp)                	# [0] return_address
	ldl	$9,376($sp)               	# [1] _temp_gra_spill6
#	.body
#	.label_state 1
	ldl	$10,344($sp)              	# [1] _temp_gra_spill2
	ldi	$gp,0($gp)                	!gpdisp!12	# [1] 0
#	.body
	ldl	$11,496($sp)              	# [2] _temp_gra_spill21
	ldl	$12,472($sp)              	# [2] _temp_gra_spill18
	ldl	$13,448($sp)              	# [3] _temp_gra_spill15
	ldl	$14,424($sp)              	# [3] _temp_gra_spill12
	ldl	$15,360($sp)              	# [4] _temp_gra_spill4
	ldi	$sp,736($sp)              	# [4] 
	ret	$31,($26),1               	# [4] 
	.align	4
.BB76_fd4t10s_sw:
#<loop> Part of loop body line 65536, head labeled .Lt_0_23554
	.loc	1	234	0
	ldl	$12,616($sp)              	# [0] _temp_gra_spill36
#	.restore $sp
	cmpeq	$27,$12,$12             	# [3] 
	.loc	1	267	0
	br	$31,.Lt_0_23810            	# [3] 
.L_CC_fd4t10s_sw:
# PU cycle count: 0.000000
	.end	fd4t10s_sw

	.section .rodata
	.org 0x0
	.align	0
	.globl	d
	.type	d, @object
	.size	d, 4
d:	# 0x0
	# offset 0
	.long	6
	# end of initialization for d
	.org 0x4
	.align	0
	# offset 4
	.ascii	"\x5f\x5a\xc4\x3f"	# float 1.53401
	.org 0x8
	.align	0
	# offset 8
	.ascii	"\x6d\xf0\xe4\x3f"	# float 1.78859
	.org 0xc
	.align	0
	# offset 12
	.ascii	"\x63\x1a\xa2\xbe"	# float -0.316608
	.org 0x10
	.align	0
	# offset 16
	.ascii	"\xb6\xe5\x9b\x3d"	# float 0.0761217
	.org 0x14
	.align	0
	# offset 20
	.ascii	"\x93\x34\x85\xbc"	# float -0.0162604
	.org 0x18
	.align	0
	# offset 24
	.ascii	"\x44\x0a\x0e\x3b"	# float 0.00216736
	.org 0x20
	.align	0
	# offset 32
	.ascii	"\x00\x00\x00\x00\x00\x00\x10\xc0"	# double -4.00000
	.org 0x28
	.align	0
	# offset 40
	.ascii	"\x00\x00\x80\xc0"	# float -4.00000
	.org 0x2c
	.align	0
	# offset 44
	.ascii	"\xab\xaa\xaa\x3d"	# float 0.0833333
	.org 0x30
	.align	0
	# offset 48
	.ascii	"\x00\x00\x00\x00\x00\x00\x00\x40"	# double 2.00000
	.weak	_tdata_local_start
	.weak	_tdata_local_end
	.weak	_tdata_private_start
	.weak	_tdata_private_end
	.weak	_tdata_local_fix_end
	.section .tdata
	.align	0
	.section .tdata_local_fix
	.align	5
	.section .text1
	.align	4
	.section .rodata
	.align	3
#	.gpvalue 0

	.section .debug_info, "", "progbits"
	.align	0
	.byte	0x67, 0x00, 0x00, 0x00, 0x02, 0x00
	.long	.debug_abbrev
	.long	0x64660108, 0x30317434, 0x6c732d73, 0x2e657661
	.long	0x706f0063, 0x43436e65, 0x342e3520, 0x732d3132
	.byte	0x77, 0x2d, 0x34, 0x38, 0x35, 0x00, 0x01, 0x00
	.long	.debug_line
	.long	0x660c0102, 0x31743464, 0x735f7330, 0x01010077
	.byte	0x04, 0x92, 0x1e, 0xe0, 0x05
	.quad	.BB1_fd4t10s_sw
	.quad	.L_CC_fd4t10s_sw
	.long	0x700c0103, 0x6d617261, 0xa0910300, 0x0000007a

	.section .debug_frame, "", "progbits"
	.align	0

	.section .debug_aranges, "", "progbits"
	.align	0
	.byte	0x2c, 0x00, 0x00, 0x00, 0x02, 0x00
	.long	.debug_info
	.byte	0x08, 0x00, 0x00, 0x00, 0x00, 0x00
	.quad	.BB1_fd4t10s_sw
	.quad	.L_CC_fd4t10s_sw - .BB1_fd4t10s_sw
	.long	0x00000000, 0x00000000, 0x00000000, 0x00000000

	.section .debug_pubnames, "", "progbits"
	.align	0
	.byte	0x1d, 0x00, 0x00, 0x00, 0x02, 0x00
	.long	.debug_info
	.long	0x0000006b, 0x00000036, 0x74346466, 0x5f733031
	.byte	0x73, 0x77, 0x00, 0x00, 0x00, 0x00, 0x00

	.section .eh_frame, "a", "progbits"
	.align	0

.LEHCIE:
	.long	.LEHCIE_end - .LEHCIE_begin
.LEHCIE_begin:
	.long 0x0
	.byte	0x01, 0x00, 0x01, 0x78, 0x1a, 0x0c, 0x1e, 0x00
	.align 3
.LEHCIE_end:

	.section .debug_abbrev, "", "progbits"
	.align	0
	.long	0x03011101, 0x13082508, 0x100b420b, 0x02000006
	.long	0x0b3a012e, 0x08030b3b, 0x408b0c3f, 0x110a400c
	.long	0x00011201, 0x00050300, 0x0b3b0b3a, 0x0a020803
	.byte	0x00, 0x00, 0x00, 0x00
	.section	.note.GNU-stack,"",@progbits
	.ident	"#SWCC Version 5.421-sw-485 : src/modeling/fd4t10s-slave.c compiled with : -O2 -OPT:IEEE_arith=2 -OPT:IEEE_arith=2 -msimd "

