; ModuleID = './linear-algebra/blas/gesummv/gesummv.c'
source_filename = "./linear-algebra/blas/gesummv/gesummv.c"
target datalayout = "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-f80:128-n8:16:32:64-S128"
target triple = "x86_64-unknown-linux-gnu"

@stderr = external local_unnamed_addr global ptr, align 8
@.str = private unnamed_addr constant [23 x i8] c"==BEGIN DUMP_ARRAYS==\0A\00", align 1
@.str.1 = private unnamed_addr constant [15 x i8] c"begin dump: %s\00", align 1
@.str.2 = private unnamed_addr constant [2 x i8] c"y\00", align 1
@.str.4 = private unnamed_addr constant [8 x i8] c"%0.2lf \00", align 1
@.str.5 = private unnamed_addr constant [17 x i8] c"\0Aend   dump: %s\0A\00", align 1
@.str.6 = private unnamed_addr constant [23 x i8] c"==END   DUMP_ARRAYS==\0A\00", align 1

; Function Attrs: nounwind uwtable
define dso_local noundef i32 @main(i32 noundef %0, ptr nocapture noundef readnone %1) local_unnamed_addr #0 {
  %3 = alloca { ptr, ptr, ptr }, align 8
  %4 = alloca { ptr, ptr, ptr }, align 8
  %5 = tail call ptr @polybench_alloc_data(i64 noundef 400000000, i32 noundef 8) #9
  %6 = ptrtoint ptr %5 to i64
  %7 = tail call ptr @polybench_alloc_data(i64 noundef 400000000, i32 noundef 8) #9
  %8 = ptrtoint ptr %7 to i64
  %9 = tail call ptr @polybench_alloc_data(i64 noundef 20000, i32 noundef 8) #9
  %10 = tail call ptr @polybench_alloc_data(i64 noundef 20000, i32 noundef 8) #9
  %11 = tail call ptr @polybench_alloc_data(i64 noundef 20000, i32 noundef 8) #9
  %12 = sub i64 %8, %6
  %13 = icmp ult i64 %12, 16
  br label %14

14:                                               ; preds = %57, %2
  %15 = phi i64 [ 0, %2 ], [ %58, %57 ]
  %16 = trunc i64 %15 to i32
  %17 = sitofp i32 %16 to double
  %18 = fdiv double %17, 2.000000e+04
  %19 = getelementptr inbounds double, ptr %10, i64 %15
  store double %18, ptr %19, align 8, !tbaa !6
  br i1 %13, label %41, label %20

20:                                               ; preds = %14
  %21 = insertelement <2 x i64> poison, i64 %15, i64 0
  %22 = shufflevector <2 x i64> %21, <2 x i64> poison, <2 x i32> zeroinitializer
  br label %23

23:                                               ; preds = %23, %20
  %24 = phi i64 [ 0, %20 ], [ %38, %23 ]
  %25 = phi <2 x i64> [ <i64 0, i64 1>, %20 ], [ %39, %23 ]
  %26 = mul nuw nsw <2 x i64> %25, %22
  %27 = trunc <2 x i64> %26 to <2 x i32>
  %28 = add <2 x i32> %27, <i32 1, i32 1>
  %29 = urem <2 x i32> %28, <i32 20000, i32 20000>
  %30 = sitofp <2 x i32> %29 to <2 x double>
  %31 = fdiv <2 x double> %30, <double 2.000000e+04, double 2.000000e+04>
  %32 = getelementptr inbounds [20000 x double], ptr %5, i64 %15, i64 %24
  store <2 x double> %31, ptr %32, align 8, !tbaa !6
  %33 = add <2 x i32> %27, <i32 2, i32 2>
  %34 = urem <2 x i32> %33, <i32 20000, i32 20000>
  %35 = sitofp <2 x i32> %34 to <2 x double>
  %36 = fdiv <2 x double> %35, <double 2.000000e+04, double 2.000000e+04>
  %37 = getelementptr inbounds [20000 x double], ptr %7, i64 %15, i64 %24
  store <2 x double> %36, ptr %37, align 8, !tbaa !6
  %38 = add nuw i64 %24, 2
  %39 = add <2 x i64> %25, <i64 2, i64 2>
  %40 = icmp eq i64 %38, 20000
  br i1 %40, label %57, label %23, !llvm.loop !10

41:                                               ; preds = %14, %41
  %42 = phi i64 [ %55, %41 ], [ 0, %14 ]
  %43 = mul nuw nsw i64 %42, %15
  %44 = trunc i64 %43 to i32
  %45 = add i32 %44, 1
  %46 = urem i32 %45, 20000
  %47 = sitofp i32 %46 to double
  %48 = fdiv double %47, 2.000000e+04
  %49 = getelementptr inbounds [20000 x double], ptr %5, i64 %15, i64 %42
  store double %48, ptr %49, align 8, !tbaa !6
  %50 = add i32 %44, 2
  %51 = urem i32 %50, 20000
  %52 = sitofp i32 %51 to double
  %53 = fdiv double %52, 2.000000e+04
  %54 = getelementptr inbounds [20000 x double], ptr %7, i64 %15, i64 %42
  store double %53, ptr %54, align 8, !tbaa !6
  %55 = add nuw nsw i64 %42, 1
  %56 = icmp eq i64 %55, 20000
  br i1 %56, label %57, label %41, !llvm.loop !13

57:                                               ; preds = %23, %41
  %58 = add nuw nsw i64 %15, 1
  %59 = icmp eq i64 %58, 20000
  br i1 %59, label %60, label %14

60:                                               ; preds = %57
  tail call void (...) @polybench_timer_start() #9
  %61 = getelementptr double, ptr %9, i64 20000
  %62 = icmp ule ptr %61, %11
  %63 = getelementptr double, ptr %11, i64 20000
  %64 = icmp ule ptr %63, %9
  %65 = or i1 %62, %64
  %66 = getelementptr double, ptr %7, i64 400000000
  %67 = icmp ule ptr %66, %11
  %68 = icmp ule ptr %63, %7
  %69 = or i1 %67, %68
  %70 = and i1 %65, %69
  %71 = getelementptr double, ptr %10, i64 20000
  %72 = icmp ule ptr %71, %11
  %73 = icmp ule ptr %63, %10
  %74 = or i1 %72, %73
  %75 = and i1 %74, %70
  %76 = getelementptr double, ptr %5, i64 400000000
  %77 = icmp ule ptr %76, %11
  %78 = icmp ule ptr %63, %5
  %79 = or i1 %77, %78
  %80 = and i1 %79, %75
  %81 = icmp ule ptr %66, %9
  %82 = icmp ule ptr %61, %7
  %83 = or i1 %81, %82
  %84 = and i1 %83, %80
  %85 = icmp ule ptr %71, %9
  %86 = icmp ule ptr %61, %10
  %87 = or i1 %86, %85
  %88 = and i1 %87, %84
  %89 = icmp ule ptr %76, %9
  %90 = icmp ule ptr %61, %5
  %91 = or i1 %89, %90
  %92 = and i1 %91, %88
  br i1 %92, label %93, label %117

93:                                               ; preds = %60
  tail call void @llvm.memset.p0.i64(ptr noundef nonnull align 8 dereferenceable(160000) %11, i8 0, i64 160000, i1 false), !alias.scope !14, !noalias !17
  store ptr %7, ptr %4, align 8
  %94 = getelementptr inbounds { ptr, ptr, ptr }, ptr %4, i64 0, i32 1
  store ptr %10, ptr %94, align 8
  %95 = getelementptr inbounds { ptr, ptr, ptr }, ptr %4, i64 0, i32 2
  store ptr %11, ptr %95, align 8
  call void @GOMP_parallel_loop_runtime_start(ptr nonnull @main_polly_subfn, ptr nonnull %4, i32 0, i64 0, i64 625, i64 1) #9
  call void @main_polly_subfn(ptr nonnull %4) #9
  call void @GOMP_parallel_end() #9
  call void @llvm.memset.p0.i64(ptr noundef nonnull align 8 dereferenceable(160000) %9, i8 0, i64 160000, i1 false), !alias.scope !22, !noalias !23
  store ptr %5, ptr %3, align 8
  %96 = getelementptr inbounds { ptr, ptr, ptr }, ptr %3, i64 0, i32 1
  store ptr %10, ptr %96, align 8
  %97 = getelementptr inbounds { ptr, ptr, ptr }, ptr %3, i64 0, i32 2
  store ptr %9, ptr %97, align 8
  call void @GOMP_parallel_loop_runtime_start(ptr nonnull @main_polly_subfn_7, ptr nonnull %3, i32 0, i64 0, i64 625, i64 1) #9
  call void @main_polly_subfn_7(ptr nonnull %3) #9
  call void @GOMP_parallel_end() #9
  br label %98

98:                                               ; preds = %98, %93
  %99 = phi i64 [ 0, %93 ], [ %115, %98 ]
  %100 = getelementptr double, ptr %11, i64 %99
  %101 = getelementptr double, ptr %100, i64 2
  %102 = load <2 x double>, ptr %100, align 8, !alias.scope !14, !noalias !17, !llvm.access.group !24
  %103 = load <2 x double>, ptr %101, align 8, !alias.scope !14, !noalias !17, !llvm.access.group !24
  %104 = shl nuw nsw i64 %99, 3
  %105 = getelementptr i8, ptr %9, i64 %104
  %106 = getelementptr double, ptr %105, i64 2
  %107 = load <2 x double>, ptr %105, align 8, !alias.scope !22, !noalias !23, !llvm.access.group !24
  %108 = load <2 x double>, ptr %106, align 8, !alias.scope !22, !noalias !23, !llvm.access.group !24
  %109 = fmul <2 x double> %102, <double 1.200000e+00, double 1.200000e+00>
  %110 = fmul <2 x double> %103, <double 1.200000e+00, double 1.200000e+00>
  %111 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %107, <2 x double> <double 1.500000e+00, double 1.500000e+00>, <2 x double> %109)
  %112 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %108, <2 x double> <double 1.500000e+00, double 1.500000e+00>, <2 x double> %110)
  %113 = getelementptr i8, ptr %11, i64 %104
  %114 = getelementptr double, ptr %113, i64 2
  store <2 x double> %111, ptr %113, align 8, !alias.scope !14, !noalias !17, !llvm.access.group !24
  store <2 x double> %112, ptr %114, align 8, !alias.scope !14, !noalias !17, !llvm.access.group !24
  %115 = add nuw i64 %99, 4
  %116 = icmp eq i64 %115, 20000
  br i1 %116, label %142, label %98, !llvm.loop !25

117:                                              ; preds = %60, %136
  %118 = phi i64 [ %140, %136 ], [ 0, %60 ]
  %119 = getelementptr inbounds double, ptr %9, i64 %118
  store double 0.000000e+00, ptr %119, align 8, !tbaa !6
  %120 = getelementptr inbounds double, ptr %11, i64 %118
  store double 0.000000e+00, ptr %120, align 8, !tbaa !6
  br label %121

121:                                              ; preds = %121, %117
  %122 = phi i64 [ 0, %117 ], [ %134, %121 ]
  %123 = getelementptr inbounds [20000 x double], ptr %5, i64 %118, i64 %122
  %124 = load double, ptr %123, align 8, !tbaa !6
  %125 = getelementptr inbounds double, ptr %10, i64 %122
  %126 = load double, ptr %125, align 8, !tbaa !6
  %127 = load double, ptr %119, align 8, !tbaa !6
  %128 = tail call double @llvm.fmuladd.f64(double %124, double %126, double %127)
  store double %128, ptr %119, align 8, !tbaa !6
  %129 = getelementptr inbounds [20000 x double], ptr %7, i64 %118, i64 %122
  %130 = load double, ptr %129, align 8, !tbaa !6
  %131 = load double, ptr %125, align 8, !tbaa !6
  %132 = load double, ptr %120, align 8, !tbaa !6
  %133 = tail call double @llvm.fmuladd.f64(double %130, double %131, double %132)
  store double %133, ptr %120, align 8, !tbaa !6
  %134 = add nuw nsw i64 %122, 1
  %135 = icmp eq i64 %134, 20000
  br i1 %135, label %136, label %121

136:                                              ; preds = %121
  %137 = load double, ptr %119, align 8, !tbaa !6
  %138 = fmul double %133, 1.200000e+00
  %139 = tail call double @llvm.fmuladd.f64(double %137, double 1.500000e+00, double %138)
  store double %139, ptr %120, align 8, !tbaa !6
  %140 = add nuw nsw i64 %118, 1
  %141 = icmp eq i64 %140, 20000
  br i1 %141, label %142, label %117

142:                                              ; preds = %136, %98
  tail call void (...) @polybench_timer_stop() #9
  tail call void (...) @polybench_timer_print() #9
  %143 = load ptr, ptr @stderr, align 8, !tbaa !27
  %144 = tail call i64 @fwrite(ptr nonnull @.str, i64 22, i64 1, ptr %143) #10
  %145 = load ptr, ptr @stderr, align 8, !tbaa !27
  %146 = tail call i32 (ptr, ptr, ...) @fprintf(ptr noundef %145, ptr noundef nonnull @.str.1, ptr noundef nonnull @.str.2) #10
  br label %147

147:                                              ; preds = %155, %142
  %148 = phi i64 [ 0, %142 ], [ %160, %155 ]
  %149 = trunc i64 %148 to i16
  %150 = urem i16 %149, 20
  %151 = icmp eq i16 %150, 0
  br i1 %151, label %152, label %155

152:                                              ; preds = %147
  %153 = load ptr, ptr @stderr, align 8, !tbaa !27
  %154 = tail call i32 @fputc(i32 10, ptr %153)
  br label %155

155:                                              ; preds = %152, %147
  %156 = load ptr, ptr @stderr, align 8, !tbaa !27
  %157 = getelementptr inbounds double, ptr %11, i64 %148
  %158 = load double, ptr %157, align 8, !tbaa !6
  %159 = tail call i32 (ptr, ptr, ...) @fprintf(ptr noundef %156, ptr noundef nonnull @.str.4, double noundef %158) #10
  %160 = add nuw nsw i64 %148, 1
  %161 = icmp eq i64 %160, 20000
  br i1 %161, label %162, label %147

162:                                              ; preds = %155
  %163 = load ptr, ptr @stderr, align 8, !tbaa !27
  %164 = tail call i32 (ptr, ptr, ...) @fprintf(ptr noundef %163, ptr noundef nonnull @.str.5, ptr noundef nonnull @.str.2) #10
  %165 = load ptr, ptr @stderr, align 8, !tbaa !27
  %166 = tail call i64 @fwrite(ptr nonnull @.str.6, i64 22, i64 1, ptr %165) #10
  tail call void @free(ptr noundef %5) #9
  tail call void @free(ptr noundef %7) #9
  tail call void @free(ptr noundef %9) #9
  tail call void @free(ptr noundef %10) #9
  tail call void @free(ptr noundef nonnull %11) #9
  ret i32 0
}

declare ptr @polybench_alloc_data(i64 noundef, i32 noundef) local_unnamed_addr #1

declare void @polybench_timer_start(...) local_unnamed_addr #1

declare void @polybench_timer_stop(...) local_unnamed_addr #1

declare void @polybench_timer_print(...) local_unnamed_addr #1

; Function Attrs: mustprogress nounwind willreturn allockind("free") memory(argmem: readwrite, inaccessiblemem: readwrite)
declare void @free(ptr allocptr nocapture noundef) local_unnamed_addr #2

; Function Attrs: mustprogress nocallback nofree nosync nounwind speculatable willreturn memory(none)
declare double @llvm.fmuladd.f64(double, double, double) #3

; Function Attrs: nofree nounwind
declare noundef i32 @fprintf(ptr nocapture noundef, ptr nocapture noundef readonly, ...) local_unnamed_addr #4

; Function Attrs: nofree nounwind
declare noundef i64 @fwrite(ptr nocapture noundef, i64 noundef, i64 noundef, ptr nocapture noundef) local_unnamed_addr #5

; Function Attrs: nofree nounwind
declare noundef i32 @fputc(i32 noundef, ptr nocapture noundef) local_unnamed_addr #5

define internal void @main_polly_subfn(ptr %0) #6 {
  %2 = alloca i64, align 8
  %3 = alloca i64, align 8
  %4 = load ptr, ptr %0, align 8
  %5 = getelementptr inbounds { ptr, ptr, ptr }, ptr %0, i64 0, i32 1
  %6 = load ptr, ptr %5, align 8
  %7 = getelementptr inbounds { ptr, ptr, ptr }, ptr %0, i64 0, i32 2
  %8 = load ptr, ptr %7, align 8
  %9 = call i8 @GOMP_loop_runtime_next(ptr nonnull %2, ptr nonnull %3)
  %10 = icmp eq i8 %9, 0
  br i1 %10, label %11, label %18

11:                                               ; preds = %12, %1
  call void @GOMP_loop_end_nowait()
  ret void

12:                                               ; preds = %15
  %13 = call i8 @GOMP_loop_runtime_next(ptr nonnull %2, ptr nonnull %3)
  %14 = icmp eq i8 %13, 0
  br i1 %14, label %11, label %18

15:                                               ; preds = %425
  %16 = add nsw i64 %23, 1
  %17 = icmp slt i64 %23, %21
  br i1 %17, label %22, label %12

18:                                               ; preds = %1, %12
  %19 = load i64, ptr %2, align 8
  %20 = load i64, ptr %3, align 8
  %21 = add i64 %20, -1
  br label %22

22:                                               ; preds = %18, %15
  %23 = phi i64 [ %19, %18 ], [ %16, %15 ]
  %24 = shl nsw i64 %23, 5
  br label %25

25:                                               ; preds = %22, %425
  %26 = phi i64 [ 0, %22 ], [ %426, %425 ]
  %27 = shl i64 %26, 8
  %28 = or disjoint i64 %27, 8
  %29 = or disjoint i64 %27, 16
  %30 = or disjoint i64 %27, 24
  %31 = or disjoint i64 %27, 32
  %32 = or disjoint i64 %27, 40
  %33 = or disjoint i64 %27, 48
  %34 = or disjoint i64 %27, 56
  %35 = or disjoint i64 %27, 64
  %36 = or disjoint i64 %27, 72
  %37 = or disjoint i64 %27, 80
  %38 = or disjoint i64 %27, 88
  %39 = or disjoint i64 %27, 96
  %40 = or disjoint i64 %27, 104
  %41 = or disjoint i64 %27, 112
  %42 = or disjoint i64 %27, 120
  %43 = or disjoint i64 %27, 128
  %44 = or disjoint i64 %27, 136
  %45 = or disjoint i64 %27, 144
  %46 = or disjoint i64 %27, 152
  %47 = or disjoint i64 %27, 160
  %48 = or disjoint i64 %27, 168
  %49 = or disjoint i64 %27, 176
  %50 = or disjoint i64 %27, 184
  %51 = or disjoint i64 %27, 192
  %52 = or disjoint i64 %27, 200
  %53 = or disjoint i64 %27, 208
  %54 = or disjoint i64 %27, 216
  %55 = or disjoint i64 %27, 224
  %56 = or disjoint i64 %27, 232
  %57 = or disjoint i64 %27, 240
  %58 = or disjoint i64 %27, 248
  %59 = getelementptr i8, ptr %6, i64 %58
  %60 = load double, ptr %59, align 8, !alias.scope !29, !noalias !30
  %61 = getelementptr i8, ptr %6, i64 %57
  %62 = load double, ptr %61, align 8, !alias.scope !29, !noalias !30
  %63 = getelementptr i8, ptr %6, i64 %56
  %64 = load double, ptr %63, align 8, !alias.scope !29, !noalias !30
  %65 = getelementptr i8, ptr %6, i64 %55
  %66 = load double, ptr %65, align 8, !alias.scope !29, !noalias !30
  %67 = getelementptr i8, ptr %6, i64 %54
  %68 = load double, ptr %67, align 8, !alias.scope !29, !noalias !30
  %69 = getelementptr i8, ptr %6, i64 %53
  %70 = load double, ptr %69, align 8, !alias.scope !29, !noalias !30
  %71 = getelementptr i8, ptr %6, i64 %52
  %72 = load double, ptr %71, align 8, !alias.scope !29, !noalias !30
  %73 = getelementptr i8, ptr %6, i64 %51
  %74 = load double, ptr %73, align 8, !alias.scope !29, !noalias !30
  %75 = getelementptr i8, ptr %6, i64 %50
  %76 = load double, ptr %75, align 8, !alias.scope !29, !noalias !30
  %77 = getelementptr i8, ptr %6, i64 %49
  %78 = load double, ptr %77, align 8, !alias.scope !29, !noalias !30
  %79 = getelementptr i8, ptr %6, i64 %48
  %80 = load double, ptr %79, align 8, !alias.scope !29, !noalias !30
  %81 = getelementptr i8, ptr %6, i64 %47
  %82 = load double, ptr %81, align 8, !alias.scope !29, !noalias !30
  %83 = getelementptr i8, ptr %6, i64 %46
  %84 = load double, ptr %83, align 8, !alias.scope !29, !noalias !30
  %85 = getelementptr i8, ptr %6, i64 %45
  %86 = load double, ptr %85, align 8, !alias.scope !29, !noalias !30
  %87 = getelementptr i8, ptr %6, i64 %44
  %88 = load double, ptr %87, align 8, !alias.scope !29, !noalias !30
  %89 = getelementptr i8, ptr %6, i64 %43
  %90 = load double, ptr %89, align 8, !alias.scope !29, !noalias !30
  %91 = getelementptr i8, ptr %6, i64 %42
  %92 = load double, ptr %91, align 8, !alias.scope !29, !noalias !30
  %93 = getelementptr i8, ptr %6, i64 %41
  %94 = load double, ptr %93, align 8, !alias.scope !29, !noalias !30
  %95 = getelementptr i8, ptr %6, i64 %40
  %96 = load double, ptr %95, align 8, !alias.scope !29, !noalias !30
  %97 = getelementptr i8, ptr %6, i64 %39
  %98 = load double, ptr %97, align 8, !alias.scope !29, !noalias !30
  %99 = getelementptr i8, ptr %6, i64 %38
  %100 = load double, ptr %99, align 8, !alias.scope !29, !noalias !30
  %101 = getelementptr i8, ptr %6, i64 %37
  %102 = load double, ptr %101, align 8, !alias.scope !29, !noalias !30
  %103 = getelementptr i8, ptr %6, i64 %36
  %104 = load double, ptr %103, align 8, !alias.scope !29, !noalias !30
  %105 = getelementptr i8, ptr %6, i64 %35
  %106 = load double, ptr %105, align 8, !alias.scope !29, !noalias !30
  %107 = getelementptr i8, ptr %6, i64 %34
  %108 = load double, ptr %107, align 8, !alias.scope !29, !noalias !30
  %109 = getelementptr i8, ptr %6, i64 %33
  %110 = load double, ptr %109, align 8, !alias.scope !29, !noalias !30
  %111 = getelementptr i8, ptr %6, i64 %32
  %112 = load double, ptr %111, align 8, !alias.scope !29, !noalias !30
  %113 = getelementptr i8, ptr %6, i64 %31
  %114 = load double, ptr %113, align 8, !alias.scope !29, !noalias !30
  %115 = getelementptr i8, ptr %6, i64 %30
  %116 = load double, ptr %115, align 8, !alias.scope !29, !noalias !30
  %117 = getelementptr i8, ptr %6, i64 %29
  %118 = load double, ptr %117, align 8, !alias.scope !29, !noalias !30
  %119 = getelementptr i8, ptr %6, i64 %28
  %120 = load double, ptr %119, align 8, !alias.scope !29, !noalias !30
  %121 = getelementptr i8, ptr %6, i64 %27
  %122 = load double, ptr %121, align 8, !alias.scope !29, !noalias !30
  %123 = insertelement <2 x double> poison, double %122, i64 0
  %124 = shufflevector <2 x double> %123, <2 x double> poison, <2 x i32> zeroinitializer
  %125 = insertelement <2 x double> poison, double %120, i64 0
  %126 = shufflevector <2 x double> %125, <2 x double> poison, <2 x i32> zeroinitializer
  %127 = insertelement <2 x double> poison, double %118, i64 0
  %128 = shufflevector <2 x double> %127, <2 x double> poison, <2 x i32> zeroinitializer
  %129 = insertelement <2 x double> poison, double %116, i64 0
  %130 = shufflevector <2 x double> %129, <2 x double> poison, <2 x i32> zeroinitializer
  %131 = insertelement <2 x double> poison, double %114, i64 0
  %132 = shufflevector <2 x double> %131, <2 x double> poison, <2 x i32> zeroinitializer
  %133 = insertelement <2 x double> poison, double %112, i64 0
  %134 = shufflevector <2 x double> %133, <2 x double> poison, <2 x i32> zeroinitializer
  %135 = insertelement <2 x double> poison, double %110, i64 0
  %136 = shufflevector <2 x double> %135, <2 x double> poison, <2 x i32> zeroinitializer
  %137 = insertelement <2 x double> poison, double %108, i64 0
  %138 = shufflevector <2 x double> %137, <2 x double> poison, <2 x i32> zeroinitializer
  %139 = insertelement <2 x double> poison, double %106, i64 0
  %140 = shufflevector <2 x double> %139, <2 x double> poison, <2 x i32> zeroinitializer
  %141 = insertelement <2 x double> poison, double %104, i64 0
  %142 = shufflevector <2 x double> %141, <2 x double> poison, <2 x i32> zeroinitializer
  %143 = insertelement <2 x double> poison, double %102, i64 0
  %144 = shufflevector <2 x double> %143, <2 x double> poison, <2 x i32> zeroinitializer
  %145 = insertelement <2 x double> poison, double %100, i64 0
  %146 = shufflevector <2 x double> %145, <2 x double> poison, <2 x i32> zeroinitializer
  %147 = insertelement <2 x double> poison, double %98, i64 0
  %148 = shufflevector <2 x double> %147, <2 x double> poison, <2 x i32> zeroinitializer
  %149 = insertelement <2 x double> poison, double %96, i64 0
  %150 = shufflevector <2 x double> %149, <2 x double> poison, <2 x i32> zeroinitializer
  %151 = insertelement <2 x double> poison, double %94, i64 0
  %152 = shufflevector <2 x double> %151, <2 x double> poison, <2 x i32> zeroinitializer
  %153 = insertelement <2 x double> poison, double %92, i64 0
  %154 = shufflevector <2 x double> %153, <2 x double> poison, <2 x i32> zeroinitializer
  %155 = insertelement <2 x double> poison, double %90, i64 0
  %156 = shufflevector <2 x double> %155, <2 x double> poison, <2 x i32> zeroinitializer
  %157 = insertelement <2 x double> poison, double %88, i64 0
  %158 = shufflevector <2 x double> %157, <2 x double> poison, <2 x i32> zeroinitializer
  %159 = insertelement <2 x double> poison, double %86, i64 0
  %160 = shufflevector <2 x double> %159, <2 x double> poison, <2 x i32> zeroinitializer
  %161 = insertelement <2 x double> poison, double %84, i64 0
  %162 = shufflevector <2 x double> %161, <2 x double> poison, <2 x i32> zeroinitializer
  %163 = insertelement <2 x double> poison, double %82, i64 0
  %164 = shufflevector <2 x double> %163, <2 x double> poison, <2 x i32> zeroinitializer
  %165 = insertelement <2 x double> poison, double %80, i64 0
  %166 = shufflevector <2 x double> %165, <2 x double> poison, <2 x i32> zeroinitializer
  %167 = insertelement <2 x double> poison, double %78, i64 0
  %168 = shufflevector <2 x double> %167, <2 x double> poison, <2 x i32> zeroinitializer
  %169 = insertelement <2 x double> poison, double %76, i64 0
  %170 = shufflevector <2 x double> %169, <2 x double> poison, <2 x i32> zeroinitializer
  %171 = insertelement <2 x double> poison, double %74, i64 0
  %172 = shufflevector <2 x double> %171, <2 x double> poison, <2 x i32> zeroinitializer
  %173 = insertelement <2 x double> poison, double %72, i64 0
  %174 = shufflevector <2 x double> %173, <2 x double> poison, <2 x i32> zeroinitializer
  %175 = insertelement <2 x double> poison, double %70, i64 0
  %176 = shufflevector <2 x double> %175, <2 x double> poison, <2 x i32> zeroinitializer
  %177 = insertelement <2 x double> poison, double %68, i64 0
  %178 = shufflevector <2 x double> %177, <2 x double> poison, <2 x i32> zeroinitializer
  %179 = insertelement <2 x double> poison, double %66, i64 0
  %180 = shufflevector <2 x double> %179, <2 x double> poison, <2 x i32> zeroinitializer
  %181 = insertelement <2 x double> poison, double %64, i64 0
  %182 = shufflevector <2 x double> %181, <2 x double> poison, <2 x i32> zeroinitializer
  %183 = insertelement <2 x double> poison, double %62, i64 0
  %184 = shufflevector <2 x double> %183, <2 x double> poison, <2 x i32> zeroinitializer
  %185 = insertelement <2 x double> poison, double %60, i64 0
  %186 = shufflevector <2 x double> %185, <2 x double> poison, <2 x i32> zeroinitializer
  br label %187

187:                                              ; preds = %187, %25
  %188 = phi i64 [ 0, %25 ], [ %423, %187 ]
  %189 = or disjoint i64 %188, 1
  %190 = add nsw i64 %188, %24
  %191 = add nsw i64 %189, %24
  %192 = mul i64 %190, 160000
  %193 = mul i64 %191, 160000
  %194 = getelementptr i8, ptr %4, i64 %192
  %195 = getelementptr i8, ptr %4, i64 %193
  %196 = shl i64 %190, 3
  %197 = getelementptr i8, ptr %8, i64 %196
  %198 = load <2 x double>, ptr %197, align 8, !alias.scope !14, !noalias !17
  %199 = getelementptr i8, ptr %194, i64 %27
  %200 = getelementptr i8, ptr %195, i64 %27
  %201 = load double, ptr %199, align 8, !alias.scope !31, !noalias !32
  %202 = load double, ptr %200, align 8, !alias.scope !31, !noalias !32
  %203 = insertelement <2 x double> poison, double %201, i64 0
  %204 = insertelement <2 x double> %203, double %202, i64 1
  %205 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %204, <2 x double> %124, <2 x double> %198)
  %206 = getelementptr i8, ptr %194, i64 %28
  %207 = getelementptr i8, ptr %195, i64 %28
  %208 = load double, ptr %206, align 8, !alias.scope !31, !noalias !32
  %209 = load double, ptr %207, align 8, !alias.scope !31, !noalias !32
  %210 = insertelement <2 x double> poison, double %208, i64 0
  %211 = insertelement <2 x double> %210, double %209, i64 1
  %212 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %211, <2 x double> %126, <2 x double> %205)
  %213 = getelementptr i8, ptr %194, i64 %29
  %214 = getelementptr i8, ptr %195, i64 %29
  %215 = load double, ptr %213, align 8, !alias.scope !31, !noalias !32
  %216 = load double, ptr %214, align 8, !alias.scope !31, !noalias !32
  %217 = insertelement <2 x double> poison, double %215, i64 0
  %218 = insertelement <2 x double> %217, double %216, i64 1
  %219 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %218, <2 x double> %128, <2 x double> %212)
  %220 = getelementptr i8, ptr %194, i64 %30
  %221 = getelementptr i8, ptr %195, i64 %30
  %222 = load double, ptr %220, align 8, !alias.scope !31, !noalias !32
  %223 = load double, ptr %221, align 8, !alias.scope !31, !noalias !32
  %224 = insertelement <2 x double> poison, double %222, i64 0
  %225 = insertelement <2 x double> %224, double %223, i64 1
  %226 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %225, <2 x double> %130, <2 x double> %219)
  %227 = getelementptr i8, ptr %194, i64 %31
  %228 = getelementptr i8, ptr %195, i64 %31
  %229 = load double, ptr %227, align 8, !alias.scope !31, !noalias !32
  %230 = load double, ptr %228, align 8, !alias.scope !31, !noalias !32
  %231 = insertelement <2 x double> poison, double %229, i64 0
  %232 = insertelement <2 x double> %231, double %230, i64 1
  %233 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %232, <2 x double> %132, <2 x double> %226)
  %234 = getelementptr i8, ptr %194, i64 %32
  %235 = getelementptr i8, ptr %195, i64 %32
  %236 = load double, ptr %234, align 8, !alias.scope !31, !noalias !32
  %237 = load double, ptr %235, align 8, !alias.scope !31, !noalias !32
  %238 = insertelement <2 x double> poison, double %236, i64 0
  %239 = insertelement <2 x double> %238, double %237, i64 1
  %240 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %239, <2 x double> %134, <2 x double> %233)
  %241 = getelementptr i8, ptr %194, i64 %33
  %242 = getelementptr i8, ptr %195, i64 %33
  %243 = load double, ptr %241, align 8, !alias.scope !31, !noalias !32
  %244 = load double, ptr %242, align 8, !alias.scope !31, !noalias !32
  %245 = insertelement <2 x double> poison, double %243, i64 0
  %246 = insertelement <2 x double> %245, double %244, i64 1
  %247 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %246, <2 x double> %136, <2 x double> %240)
  %248 = getelementptr i8, ptr %194, i64 %34
  %249 = getelementptr i8, ptr %195, i64 %34
  %250 = load double, ptr %248, align 8, !alias.scope !31, !noalias !32
  %251 = load double, ptr %249, align 8, !alias.scope !31, !noalias !32
  %252 = insertelement <2 x double> poison, double %250, i64 0
  %253 = insertelement <2 x double> %252, double %251, i64 1
  %254 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %253, <2 x double> %138, <2 x double> %247)
  %255 = getelementptr i8, ptr %194, i64 %35
  %256 = getelementptr i8, ptr %195, i64 %35
  %257 = load double, ptr %255, align 8, !alias.scope !31, !noalias !32
  %258 = load double, ptr %256, align 8, !alias.scope !31, !noalias !32
  %259 = insertelement <2 x double> poison, double %257, i64 0
  %260 = insertelement <2 x double> %259, double %258, i64 1
  %261 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %260, <2 x double> %140, <2 x double> %254)
  %262 = getelementptr i8, ptr %194, i64 %36
  %263 = getelementptr i8, ptr %195, i64 %36
  %264 = load double, ptr %262, align 8, !alias.scope !31, !noalias !32
  %265 = load double, ptr %263, align 8, !alias.scope !31, !noalias !32
  %266 = insertelement <2 x double> poison, double %264, i64 0
  %267 = insertelement <2 x double> %266, double %265, i64 1
  %268 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %267, <2 x double> %142, <2 x double> %261)
  %269 = getelementptr i8, ptr %194, i64 %37
  %270 = getelementptr i8, ptr %195, i64 %37
  %271 = load double, ptr %269, align 8, !alias.scope !31, !noalias !32
  %272 = load double, ptr %270, align 8, !alias.scope !31, !noalias !32
  %273 = insertelement <2 x double> poison, double %271, i64 0
  %274 = insertelement <2 x double> %273, double %272, i64 1
  %275 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %274, <2 x double> %144, <2 x double> %268)
  %276 = getelementptr i8, ptr %194, i64 %38
  %277 = getelementptr i8, ptr %195, i64 %38
  %278 = load double, ptr %276, align 8, !alias.scope !31, !noalias !32
  %279 = load double, ptr %277, align 8, !alias.scope !31, !noalias !32
  %280 = insertelement <2 x double> poison, double %278, i64 0
  %281 = insertelement <2 x double> %280, double %279, i64 1
  %282 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %281, <2 x double> %146, <2 x double> %275)
  %283 = getelementptr i8, ptr %194, i64 %39
  %284 = getelementptr i8, ptr %195, i64 %39
  %285 = load double, ptr %283, align 8, !alias.scope !31, !noalias !32
  %286 = load double, ptr %284, align 8, !alias.scope !31, !noalias !32
  %287 = insertelement <2 x double> poison, double %285, i64 0
  %288 = insertelement <2 x double> %287, double %286, i64 1
  %289 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %288, <2 x double> %148, <2 x double> %282)
  %290 = getelementptr i8, ptr %194, i64 %40
  %291 = getelementptr i8, ptr %195, i64 %40
  %292 = load double, ptr %290, align 8, !alias.scope !31, !noalias !32
  %293 = load double, ptr %291, align 8, !alias.scope !31, !noalias !32
  %294 = insertelement <2 x double> poison, double %292, i64 0
  %295 = insertelement <2 x double> %294, double %293, i64 1
  %296 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %295, <2 x double> %150, <2 x double> %289)
  %297 = getelementptr i8, ptr %194, i64 %41
  %298 = getelementptr i8, ptr %195, i64 %41
  %299 = load double, ptr %297, align 8, !alias.scope !31, !noalias !32
  %300 = load double, ptr %298, align 8, !alias.scope !31, !noalias !32
  %301 = insertelement <2 x double> poison, double %299, i64 0
  %302 = insertelement <2 x double> %301, double %300, i64 1
  %303 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %302, <2 x double> %152, <2 x double> %296)
  %304 = getelementptr i8, ptr %194, i64 %42
  %305 = getelementptr i8, ptr %195, i64 %42
  %306 = load double, ptr %304, align 8, !alias.scope !31, !noalias !32
  %307 = load double, ptr %305, align 8, !alias.scope !31, !noalias !32
  %308 = insertelement <2 x double> poison, double %306, i64 0
  %309 = insertelement <2 x double> %308, double %307, i64 1
  %310 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %309, <2 x double> %154, <2 x double> %303)
  %311 = getelementptr i8, ptr %194, i64 %43
  %312 = getelementptr i8, ptr %195, i64 %43
  %313 = load double, ptr %311, align 8, !alias.scope !31, !noalias !32
  %314 = load double, ptr %312, align 8, !alias.scope !31, !noalias !32
  %315 = insertelement <2 x double> poison, double %313, i64 0
  %316 = insertelement <2 x double> %315, double %314, i64 1
  %317 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %316, <2 x double> %156, <2 x double> %310)
  %318 = getelementptr i8, ptr %194, i64 %44
  %319 = getelementptr i8, ptr %195, i64 %44
  %320 = load double, ptr %318, align 8, !alias.scope !31, !noalias !32
  %321 = load double, ptr %319, align 8, !alias.scope !31, !noalias !32
  %322 = insertelement <2 x double> poison, double %320, i64 0
  %323 = insertelement <2 x double> %322, double %321, i64 1
  %324 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %323, <2 x double> %158, <2 x double> %317)
  %325 = getelementptr i8, ptr %194, i64 %45
  %326 = getelementptr i8, ptr %195, i64 %45
  %327 = load double, ptr %325, align 8, !alias.scope !31, !noalias !32
  %328 = load double, ptr %326, align 8, !alias.scope !31, !noalias !32
  %329 = insertelement <2 x double> poison, double %327, i64 0
  %330 = insertelement <2 x double> %329, double %328, i64 1
  %331 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %330, <2 x double> %160, <2 x double> %324)
  %332 = getelementptr i8, ptr %194, i64 %46
  %333 = getelementptr i8, ptr %195, i64 %46
  %334 = load double, ptr %332, align 8, !alias.scope !31, !noalias !32
  %335 = load double, ptr %333, align 8, !alias.scope !31, !noalias !32
  %336 = insertelement <2 x double> poison, double %334, i64 0
  %337 = insertelement <2 x double> %336, double %335, i64 1
  %338 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %337, <2 x double> %162, <2 x double> %331)
  %339 = getelementptr i8, ptr %194, i64 %47
  %340 = getelementptr i8, ptr %195, i64 %47
  %341 = load double, ptr %339, align 8, !alias.scope !31, !noalias !32
  %342 = load double, ptr %340, align 8, !alias.scope !31, !noalias !32
  %343 = insertelement <2 x double> poison, double %341, i64 0
  %344 = insertelement <2 x double> %343, double %342, i64 1
  %345 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %344, <2 x double> %164, <2 x double> %338)
  %346 = getelementptr i8, ptr %194, i64 %48
  %347 = getelementptr i8, ptr %195, i64 %48
  %348 = load double, ptr %346, align 8, !alias.scope !31, !noalias !32
  %349 = load double, ptr %347, align 8, !alias.scope !31, !noalias !32
  %350 = insertelement <2 x double> poison, double %348, i64 0
  %351 = insertelement <2 x double> %350, double %349, i64 1
  %352 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %351, <2 x double> %166, <2 x double> %345)
  %353 = getelementptr i8, ptr %194, i64 %49
  %354 = getelementptr i8, ptr %195, i64 %49
  %355 = load double, ptr %353, align 8, !alias.scope !31, !noalias !32
  %356 = load double, ptr %354, align 8, !alias.scope !31, !noalias !32
  %357 = insertelement <2 x double> poison, double %355, i64 0
  %358 = insertelement <2 x double> %357, double %356, i64 1
  %359 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %358, <2 x double> %168, <2 x double> %352)
  %360 = getelementptr i8, ptr %194, i64 %50
  %361 = getelementptr i8, ptr %195, i64 %50
  %362 = load double, ptr %360, align 8, !alias.scope !31, !noalias !32
  %363 = load double, ptr %361, align 8, !alias.scope !31, !noalias !32
  %364 = insertelement <2 x double> poison, double %362, i64 0
  %365 = insertelement <2 x double> %364, double %363, i64 1
  %366 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %365, <2 x double> %170, <2 x double> %359)
  %367 = getelementptr i8, ptr %194, i64 %51
  %368 = getelementptr i8, ptr %195, i64 %51
  %369 = load double, ptr %367, align 8, !alias.scope !31, !noalias !32
  %370 = load double, ptr %368, align 8, !alias.scope !31, !noalias !32
  %371 = insertelement <2 x double> poison, double %369, i64 0
  %372 = insertelement <2 x double> %371, double %370, i64 1
  %373 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %372, <2 x double> %172, <2 x double> %366)
  %374 = getelementptr i8, ptr %194, i64 %52
  %375 = getelementptr i8, ptr %195, i64 %52
  %376 = load double, ptr %374, align 8, !alias.scope !31, !noalias !32
  %377 = load double, ptr %375, align 8, !alias.scope !31, !noalias !32
  %378 = insertelement <2 x double> poison, double %376, i64 0
  %379 = insertelement <2 x double> %378, double %377, i64 1
  %380 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %379, <2 x double> %174, <2 x double> %373)
  %381 = getelementptr i8, ptr %194, i64 %53
  %382 = getelementptr i8, ptr %195, i64 %53
  %383 = load double, ptr %381, align 8, !alias.scope !31, !noalias !32
  %384 = load double, ptr %382, align 8, !alias.scope !31, !noalias !32
  %385 = insertelement <2 x double> poison, double %383, i64 0
  %386 = insertelement <2 x double> %385, double %384, i64 1
  %387 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %386, <2 x double> %176, <2 x double> %380)
  %388 = getelementptr i8, ptr %194, i64 %54
  %389 = getelementptr i8, ptr %195, i64 %54
  %390 = load double, ptr %388, align 8, !alias.scope !31, !noalias !32
  %391 = load double, ptr %389, align 8, !alias.scope !31, !noalias !32
  %392 = insertelement <2 x double> poison, double %390, i64 0
  %393 = insertelement <2 x double> %392, double %391, i64 1
  %394 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %393, <2 x double> %178, <2 x double> %387)
  %395 = getelementptr i8, ptr %194, i64 %55
  %396 = getelementptr i8, ptr %195, i64 %55
  %397 = load double, ptr %395, align 8, !alias.scope !31, !noalias !32
  %398 = load double, ptr %396, align 8, !alias.scope !31, !noalias !32
  %399 = insertelement <2 x double> poison, double %397, i64 0
  %400 = insertelement <2 x double> %399, double %398, i64 1
  %401 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %400, <2 x double> %180, <2 x double> %394)
  %402 = getelementptr i8, ptr %194, i64 %56
  %403 = getelementptr i8, ptr %195, i64 %56
  %404 = load double, ptr %402, align 8, !alias.scope !31, !noalias !32
  %405 = load double, ptr %403, align 8, !alias.scope !31, !noalias !32
  %406 = insertelement <2 x double> poison, double %404, i64 0
  %407 = insertelement <2 x double> %406, double %405, i64 1
  %408 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %407, <2 x double> %182, <2 x double> %401)
  %409 = getelementptr i8, ptr %194, i64 %57
  %410 = getelementptr i8, ptr %195, i64 %57
  %411 = load double, ptr %409, align 8, !alias.scope !31, !noalias !32
  %412 = load double, ptr %410, align 8, !alias.scope !31, !noalias !32
  %413 = insertelement <2 x double> poison, double %411, i64 0
  %414 = insertelement <2 x double> %413, double %412, i64 1
  %415 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %414, <2 x double> %184, <2 x double> %408)
  %416 = getelementptr i8, ptr %194, i64 %58
  %417 = getelementptr i8, ptr %195, i64 %58
  %418 = load double, ptr %416, align 8, !alias.scope !31, !noalias !32
  %419 = load double, ptr %417, align 8, !alias.scope !31, !noalias !32
  %420 = insertelement <2 x double> poison, double %418, i64 0
  %421 = insertelement <2 x double> %420, double %419, i64 1
  %422 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %421, <2 x double> %186, <2 x double> %415)
  store <2 x double> %422, ptr %197, align 8, !alias.scope !14, !noalias !17
  %423 = add nuw i64 %188, 2
  %424 = icmp eq i64 %423, 32
  br i1 %424, label %425, label %187, !llvm.loop !33

425:                                              ; preds = %187
  %426 = add nuw nsw i64 %26, 1
  %427 = icmp eq i64 %426, 625
  br i1 %427, label %15, label %25
}

declare i8 @GOMP_loop_runtime_next(ptr, ptr)

declare void @GOMP_loop_end_nowait()

declare void @GOMP_parallel_loop_runtime_start(ptr, ptr, i32, i64, i64, i64)

declare void @GOMP_parallel_end()

define internal void @main_polly_subfn_7(ptr %0) #6 {
  %2 = alloca i64, align 8
  %3 = alloca i64, align 8
  %4 = load ptr, ptr %0, align 8
  %5 = getelementptr inbounds { ptr, ptr, ptr }, ptr %0, i64 0, i32 1
  %6 = load ptr, ptr %5, align 8
  %7 = getelementptr inbounds { ptr, ptr, ptr }, ptr %0, i64 0, i32 2
  %8 = load ptr, ptr %7, align 8
  %9 = call i8 @GOMP_loop_runtime_next(ptr nonnull %2, ptr nonnull %3)
  %10 = icmp eq i8 %9, 0
  br i1 %10, label %11, label %18

11:                                               ; preds = %12, %1
  call void @GOMP_loop_end_nowait()
  ret void

12:                                               ; preds = %15
  %13 = call i8 @GOMP_loop_runtime_next(ptr nonnull %2, ptr nonnull %3)
  %14 = icmp eq i8 %13, 0
  br i1 %14, label %11, label %18

15:                                               ; preds = %425
  %16 = add nsw i64 %23, 1
  %17 = icmp slt i64 %23, %21
  br i1 %17, label %22, label %12

18:                                               ; preds = %1, %12
  %19 = load i64, ptr %2, align 8
  %20 = load i64, ptr %3, align 8
  %21 = add i64 %20, -1
  br label %22

22:                                               ; preds = %18, %15
  %23 = phi i64 [ %19, %18 ], [ %16, %15 ]
  %24 = shl nsw i64 %23, 5
  br label %25

25:                                               ; preds = %22, %425
  %26 = phi i64 [ 0, %22 ], [ %426, %425 ]
  %27 = shl i64 %26, 8
  %28 = or disjoint i64 %27, 8
  %29 = or disjoint i64 %27, 16
  %30 = or disjoint i64 %27, 24
  %31 = or disjoint i64 %27, 32
  %32 = or disjoint i64 %27, 40
  %33 = or disjoint i64 %27, 48
  %34 = or disjoint i64 %27, 56
  %35 = or disjoint i64 %27, 64
  %36 = or disjoint i64 %27, 72
  %37 = or disjoint i64 %27, 80
  %38 = or disjoint i64 %27, 88
  %39 = or disjoint i64 %27, 96
  %40 = or disjoint i64 %27, 104
  %41 = or disjoint i64 %27, 112
  %42 = or disjoint i64 %27, 120
  %43 = or disjoint i64 %27, 128
  %44 = or disjoint i64 %27, 136
  %45 = or disjoint i64 %27, 144
  %46 = or disjoint i64 %27, 152
  %47 = or disjoint i64 %27, 160
  %48 = or disjoint i64 %27, 168
  %49 = or disjoint i64 %27, 176
  %50 = or disjoint i64 %27, 184
  %51 = or disjoint i64 %27, 192
  %52 = or disjoint i64 %27, 200
  %53 = or disjoint i64 %27, 208
  %54 = or disjoint i64 %27, 216
  %55 = or disjoint i64 %27, 224
  %56 = or disjoint i64 %27, 232
  %57 = or disjoint i64 %27, 240
  %58 = or disjoint i64 %27, 248
  %59 = getelementptr i8, ptr %6, i64 %58
  %60 = load double, ptr %59, align 8, !alias.scope !29, !noalias !30
  %61 = getelementptr i8, ptr %6, i64 %57
  %62 = load double, ptr %61, align 8, !alias.scope !29, !noalias !30
  %63 = getelementptr i8, ptr %6, i64 %56
  %64 = load double, ptr %63, align 8, !alias.scope !29, !noalias !30
  %65 = getelementptr i8, ptr %6, i64 %55
  %66 = load double, ptr %65, align 8, !alias.scope !29, !noalias !30
  %67 = getelementptr i8, ptr %6, i64 %54
  %68 = load double, ptr %67, align 8, !alias.scope !29, !noalias !30
  %69 = getelementptr i8, ptr %6, i64 %53
  %70 = load double, ptr %69, align 8, !alias.scope !29, !noalias !30
  %71 = getelementptr i8, ptr %6, i64 %52
  %72 = load double, ptr %71, align 8, !alias.scope !29, !noalias !30
  %73 = getelementptr i8, ptr %6, i64 %51
  %74 = load double, ptr %73, align 8, !alias.scope !29, !noalias !30
  %75 = getelementptr i8, ptr %6, i64 %50
  %76 = load double, ptr %75, align 8, !alias.scope !29, !noalias !30
  %77 = getelementptr i8, ptr %6, i64 %49
  %78 = load double, ptr %77, align 8, !alias.scope !29, !noalias !30
  %79 = getelementptr i8, ptr %6, i64 %48
  %80 = load double, ptr %79, align 8, !alias.scope !29, !noalias !30
  %81 = getelementptr i8, ptr %6, i64 %47
  %82 = load double, ptr %81, align 8, !alias.scope !29, !noalias !30
  %83 = getelementptr i8, ptr %6, i64 %46
  %84 = load double, ptr %83, align 8, !alias.scope !29, !noalias !30
  %85 = getelementptr i8, ptr %6, i64 %45
  %86 = load double, ptr %85, align 8, !alias.scope !29, !noalias !30
  %87 = getelementptr i8, ptr %6, i64 %44
  %88 = load double, ptr %87, align 8, !alias.scope !29, !noalias !30
  %89 = getelementptr i8, ptr %6, i64 %43
  %90 = load double, ptr %89, align 8, !alias.scope !29, !noalias !30
  %91 = getelementptr i8, ptr %6, i64 %42
  %92 = load double, ptr %91, align 8, !alias.scope !29, !noalias !30
  %93 = getelementptr i8, ptr %6, i64 %41
  %94 = load double, ptr %93, align 8, !alias.scope !29, !noalias !30
  %95 = getelementptr i8, ptr %6, i64 %40
  %96 = load double, ptr %95, align 8, !alias.scope !29, !noalias !30
  %97 = getelementptr i8, ptr %6, i64 %39
  %98 = load double, ptr %97, align 8, !alias.scope !29, !noalias !30
  %99 = getelementptr i8, ptr %6, i64 %38
  %100 = load double, ptr %99, align 8, !alias.scope !29, !noalias !30
  %101 = getelementptr i8, ptr %6, i64 %37
  %102 = load double, ptr %101, align 8, !alias.scope !29, !noalias !30
  %103 = getelementptr i8, ptr %6, i64 %36
  %104 = load double, ptr %103, align 8, !alias.scope !29, !noalias !30
  %105 = getelementptr i8, ptr %6, i64 %35
  %106 = load double, ptr %105, align 8, !alias.scope !29, !noalias !30
  %107 = getelementptr i8, ptr %6, i64 %34
  %108 = load double, ptr %107, align 8, !alias.scope !29, !noalias !30
  %109 = getelementptr i8, ptr %6, i64 %33
  %110 = load double, ptr %109, align 8, !alias.scope !29, !noalias !30
  %111 = getelementptr i8, ptr %6, i64 %32
  %112 = load double, ptr %111, align 8, !alias.scope !29, !noalias !30
  %113 = getelementptr i8, ptr %6, i64 %31
  %114 = load double, ptr %113, align 8, !alias.scope !29, !noalias !30
  %115 = getelementptr i8, ptr %6, i64 %30
  %116 = load double, ptr %115, align 8, !alias.scope !29, !noalias !30
  %117 = getelementptr i8, ptr %6, i64 %29
  %118 = load double, ptr %117, align 8, !alias.scope !29, !noalias !30
  %119 = getelementptr i8, ptr %6, i64 %28
  %120 = load double, ptr %119, align 8, !alias.scope !29, !noalias !30
  %121 = getelementptr i8, ptr %6, i64 %27
  %122 = load double, ptr %121, align 8, !alias.scope !29, !noalias !30
  %123 = insertelement <2 x double> poison, double %122, i64 0
  %124 = shufflevector <2 x double> %123, <2 x double> poison, <2 x i32> zeroinitializer
  %125 = insertelement <2 x double> poison, double %120, i64 0
  %126 = shufflevector <2 x double> %125, <2 x double> poison, <2 x i32> zeroinitializer
  %127 = insertelement <2 x double> poison, double %118, i64 0
  %128 = shufflevector <2 x double> %127, <2 x double> poison, <2 x i32> zeroinitializer
  %129 = insertelement <2 x double> poison, double %116, i64 0
  %130 = shufflevector <2 x double> %129, <2 x double> poison, <2 x i32> zeroinitializer
  %131 = insertelement <2 x double> poison, double %114, i64 0
  %132 = shufflevector <2 x double> %131, <2 x double> poison, <2 x i32> zeroinitializer
  %133 = insertelement <2 x double> poison, double %112, i64 0
  %134 = shufflevector <2 x double> %133, <2 x double> poison, <2 x i32> zeroinitializer
  %135 = insertelement <2 x double> poison, double %110, i64 0
  %136 = shufflevector <2 x double> %135, <2 x double> poison, <2 x i32> zeroinitializer
  %137 = insertelement <2 x double> poison, double %108, i64 0
  %138 = shufflevector <2 x double> %137, <2 x double> poison, <2 x i32> zeroinitializer
  %139 = insertelement <2 x double> poison, double %106, i64 0
  %140 = shufflevector <2 x double> %139, <2 x double> poison, <2 x i32> zeroinitializer
  %141 = insertelement <2 x double> poison, double %104, i64 0
  %142 = shufflevector <2 x double> %141, <2 x double> poison, <2 x i32> zeroinitializer
  %143 = insertelement <2 x double> poison, double %102, i64 0
  %144 = shufflevector <2 x double> %143, <2 x double> poison, <2 x i32> zeroinitializer
  %145 = insertelement <2 x double> poison, double %100, i64 0
  %146 = shufflevector <2 x double> %145, <2 x double> poison, <2 x i32> zeroinitializer
  %147 = insertelement <2 x double> poison, double %98, i64 0
  %148 = shufflevector <2 x double> %147, <2 x double> poison, <2 x i32> zeroinitializer
  %149 = insertelement <2 x double> poison, double %96, i64 0
  %150 = shufflevector <2 x double> %149, <2 x double> poison, <2 x i32> zeroinitializer
  %151 = insertelement <2 x double> poison, double %94, i64 0
  %152 = shufflevector <2 x double> %151, <2 x double> poison, <2 x i32> zeroinitializer
  %153 = insertelement <2 x double> poison, double %92, i64 0
  %154 = shufflevector <2 x double> %153, <2 x double> poison, <2 x i32> zeroinitializer
  %155 = insertelement <2 x double> poison, double %90, i64 0
  %156 = shufflevector <2 x double> %155, <2 x double> poison, <2 x i32> zeroinitializer
  %157 = insertelement <2 x double> poison, double %88, i64 0
  %158 = shufflevector <2 x double> %157, <2 x double> poison, <2 x i32> zeroinitializer
  %159 = insertelement <2 x double> poison, double %86, i64 0
  %160 = shufflevector <2 x double> %159, <2 x double> poison, <2 x i32> zeroinitializer
  %161 = insertelement <2 x double> poison, double %84, i64 0
  %162 = shufflevector <2 x double> %161, <2 x double> poison, <2 x i32> zeroinitializer
  %163 = insertelement <2 x double> poison, double %82, i64 0
  %164 = shufflevector <2 x double> %163, <2 x double> poison, <2 x i32> zeroinitializer
  %165 = insertelement <2 x double> poison, double %80, i64 0
  %166 = shufflevector <2 x double> %165, <2 x double> poison, <2 x i32> zeroinitializer
  %167 = insertelement <2 x double> poison, double %78, i64 0
  %168 = shufflevector <2 x double> %167, <2 x double> poison, <2 x i32> zeroinitializer
  %169 = insertelement <2 x double> poison, double %76, i64 0
  %170 = shufflevector <2 x double> %169, <2 x double> poison, <2 x i32> zeroinitializer
  %171 = insertelement <2 x double> poison, double %74, i64 0
  %172 = shufflevector <2 x double> %171, <2 x double> poison, <2 x i32> zeroinitializer
  %173 = insertelement <2 x double> poison, double %72, i64 0
  %174 = shufflevector <2 x double> %173, <2 x double> poison, <2 x i32> zeroinitializer
  %175 = insertelement <2 x double> poison, double %70, i64 0
  %176 = shufflevector <2 x double> %175, <2 x double> poison, <2 x i32> zeroinitializer
  %177 = insertelement <2 x double> poison, double %68, i64 0
  %178 = shufflevector <2 x double> %177, <2 x double> poison, <2 x i32> zeroinitializer
  %179 = insertelement <2 x double> poison, double %66, i64 0
  %180 = shufflevector <2 x double> %179, <2 x double> poison, <2 x i32> zeroinitializer
  %181 = insertelement <2 x double> poison, double %64, i64 0
  %182 = shufflevector <2 x double> %181, <2 x double> poison, <2 x i32> zeroinitializer
  %183 = insertelement <2 x double> poison, double %62, i64 0
  %184 = shufflevector <2 x double> %183, <2 x double> poison, <2 x i32> zeroinitializer
  %185 = insertelement <2 x double> poison, double %60, i64 0
  %186 = shufflevector <2 x double> %185, <2 x double> poison, <2 x i32> zeroinitializer
  br label %187

187:                                              ; preds = %187, %25
  %188 = phi i64 [ 0, %25 ], [ %423, %187 ]
  %189 = or disjoint i64 %188, 1
  %190 = add nsw i64 %188, %24
  %191 = add nsw i64 %189, %24
  %192 = mul i64 %190, 160000
  %193 = mul i64 %191, 160000
  %194 = getelementptr i8, ptr %4, i64 %192
  %195 = getelementptr i8, ptr %4, i64 %193
  %196 = shl i64 %190, 3
  %197 = getelementptr i8, ptr %8, i64 %196
  %198 = load <2 x double>, ptr %197, align 8, !alias.scope !22, !noalias !23
  %199 = getelementptr i8, ptr %194, i64 %27
  %200 = getelementptr i8, ptr %195, i64 %27
  %201 = load double, ptr %199, align 8, !alias.scope !34, !noalias !35
  %202 = load double, ptr %200, align 8, !alias.scope !34, !noalias !35
  %203 = insertelement <2 x double> poison, double %201, i64 0
  %204 = insertelement <2 x double> %203, double %202, i64 1
  %205 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %204, <2 x double> %124, <2 x double> %198)
  %206 = getelementptr i8, ptr %194, i64 %28
  %207 = getelementptr i8, ptr %195, i64 %28
  %208 = load double, ptr %206, align 8, !alias.scope !34, !noalias !35
  %209 = load double, ptr %207, align 8, !alias.scope !34, !noalias !35
  %210 = insertelement <2 x double> poison, double %208, i64 0
  %211 = insertelement <2 x double> %210, double %209, i64 1
  %212 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %211, <2 x double> %126, <2 x double> %205)
  %213 = getelementptr i8, ptr %194, i64 %29
  %214 = getelementptr i8, ptr %195, i64 %29
  %215 = load double, ptr %213, align 8, !alias.scope !34, !noalias !35
  %216 = load double, ptr %214, align 8, !alias.scope !34, !noalias !35
  %217 = insertelement <2 x double> poison, double %215, i64 0
  %218 = insertelement <2 x double> %217, double %216, i64 1
  %219 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %218, <2 x double> %128, <2 x double> %212)
  %220 = getelementptr i8, ptr %194, i64 %30
  %221 = getelementptr i8, ptr %195, i64 %30
  %222 = load double, ptr %220, align 8, !alias.scope !34, !noalias !35
  %223 = load double, ptr %221, align 8, !alias.scope !34, !noalias !35
  %224 = insertelement <2 x double> poison, double %222, i64 0
  %225 = insertelement <2 x double> %224, double %223, i64 1
  %226 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %225, <2 x double> %130, <2 x double> %219)
  %227 = getelementptr i8, ptr %194, i64 %31
  %228 = getelementptr i8, ptr %195, i64 %31
  %229 = load double, ptr %227, align 8, !alias.scope !34, !noalias !35
  %230 = load double, ptr %228, align 8, !alias.scope !34, !noalias !35
  %231 = insertelement <2 x double> poison, double %229, i64 0
  %232 = insertelement <2 x double> %231, double %230, i64 1
  %233 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %232, <2 x double> %132, <2 x double> %226)
  %234 = getelementptr i8, ptr %194, i64 %32
  %235 = getelementptr i8, ptr %195, i64 %32
  %236 = load double, ptr %234, align 8, !alias.scope !34, !noalias !35
  %237 = load double, ptr %235, align 8, !alias.scope !34, !noalias !35
  %238 = insertelement <2 x double> poison, double %236, i64 0
  %239 = insertelement <2 x double> %238, double %237, i64 1
  %240 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %239, <2 x double> %134, <2 x double> %233)
  %241 = getelementptr i8, ptr %194, i64 %33
  %242 = getelementptr i8, ptr %195, i64 %33
  %243 = load double, ptr %241, align 8, !alias.scope !34, !noalias !35
  %244 = load double, ptr %242, align 8, !alias.scope !34, !noalias !35
  %245 = insertelement <2 x double> poison, double %243, i64 0
  %246 = insertelement <2 x double> %245, double %244, i64 1
  %247 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %246, <2 x double> %136, <2 x double> %240)
  %248 = getelementptr i8, ptr %194, i64 %34
  %249 = getelementptr i8, ptr %195, i64 %34
  %250 = load double, ptr %248, align 8, !alias.scope !34, !noalias !35
  %251 = load double, ptr %249, align 8, !alias.scope !34, !noalias !35
  %252 = insertelement <2 x double> poison, double %250, i64 0
  %253 = insertelement <2 x double> %252, double %251, i64 1
  %254 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %253, <2 x double> %138, <2 x double> %247)
  %255 = getelementptr i8, ptr %194, i64 %35
  %256 = getelementptr i8, ptr %195, i64 %35
  %257 = load double, ptr %255, align 8, !alias.scope !34, !noalias !35
  %258 = load double, ptr %256, align 8, !alias.scope !34, !noalias !35
  %259 = insertelement <2 x double> poison, double %257, i64 0
  %260 = insertelement <2 x double> %259, double %258, i64 1
  %261 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %260, <2 x double> %140, <2 x double> %254)
  %262 = getelementptr i8, ptr %194, i64 %36
  %263 = getelementptr i8, ptr %195, i64 %36
  %264 = load double, ptr %262, align 8, !alias.scope !34, !noalias !35
  %265 = load double, ptr %263, align 8, !alias.scope !34, !noalias !35
  %266 = insertelement <2 x double> poison, double %264, i64 0
  %267 = insertelement <2 x double> %266, double %265, i64 1
  %268 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %267, <2 x double> %142, <2 x double> %261)
  %269 = getelementptr i8, ptr %194, i64 %37
  %270 = getelementptr i8, ptr %195, i64 %37
  %271 = load double, ptr %269, align 8, !alias.scope !34, !noalias !35
  %272 = load double, ptr %270, align 8, !alias.scope !34, !noalias !35
  %273 = insertelement <2 x double> poison, double %271, i64 0
  %274 = insertelement <2 x double> %273, double %272, i64 1
  %275 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %274, <2 x double> %144, <2 x double> %268)
  %276 = getelementptr i8, ptr %194, i64 %38
  %277 = getelementptr i8, ptr %195, i64 %38
  %278 = load double, ptr %276, align 8, !alias.scope !34, !noalias !35
  %279 = load double, ptr %277, align 8, !alias.scope !34, !noalias !35
  %280 = insertelement <2 x double> poison, double %278, i64 0
  %281 = insertelement <2 x double> %280, double %279, i64 1
  %282 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %281, <2 x double> %146, <2 x double> %275)
  %283 = getelementptr i8, ptr %194, i64 %39
  %284 = getelementptr i8, ptr %195, i64 %39
  %285 = load double, ptr %283, align 8, !alias.scope !34, !noalias !35
  %286 = load double, ptr %284, align 8, !alias.scope !34, !noalias !35
  %287 = insertelement <2 x double> poison, double %285, i64 0
  %288 = insertelement <2 x double> %287, double %286, i64 1
  %289 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %288, <2 x double> %148, <2 x double> %282)
  %290 = getelementptr i8, ptr %194, i64 %40
  %291 = getelementptr i8, ptr %195, i64 %40
  %292 = load double, ptr %290, align 8, !alias.scope !34, !noalias !35
  %293 = load double, ptr %291, align 8, !alias.scope !34, !noalias !35
  %294 = insertelement <2 x double> poison, double %292, i64 0
  %295 = insertelement <2 x double> %294, double %293, i64 1
  %296 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %295, <2 x double> %150, <2 x double> %289)
  %297 = getelementptr i8, ptr %194, i64 %41
  %298 = getelementptr i8, ptr %195, i64 %41
  %299 = load double, ptr %297, align 8, !alias.scope !34, !noalias !35
  %300 = load double, ptr %298, align 8, !alias.scope !34, !noalias !35
  %301 = insertelement <2 x double> poison, double %299, i64 0
  %302 = insertelement <2 x double> %301, double %300, i64 1
  %303 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %302, <2 x double> %152, <2 x double> %296)
  %304 = getelementptr i8, ptr %194, i64 %42
  %305 = getelementptr i8, ptr %195, i64 %42
  %306 = load double, ptr %304, align 8, !alias.scope !34, !noalias !35
  %307 = load double, ptr %305, align 8, !alias.scope !34, !noalias !35
  %308 = insertelement <2 x double> poison, double %306, i64 0
  %309 = insertelement <2 x double> %308, double %307, i64 1
  %310 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %309, <2 x double> %154, <2 x double> %303)
  %311 = getelementptr i8, ptr %194, i64 %43
  %312 = getelementptr i8, ptr %195, i64 %43
  %313 = load double, ptr %311, align 8, !alias.scope !34, !noalias !35
  %314 = load double, ptr %312, align 8, !alias.scope !34, !noalias !35
  %315 = insertelement <2 x double> poison, double %313, i64 0
  %316 = insertelement <2 x double> %315, double %314, i64 1
  %317 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %316, <2 x double> %156, <2 x double> %310)
  %318 = getelementptr i8, ptr %194, i64 %44
  %319 = getelementptr i8, ptr %195, i64 %44
  %320 = load double, ptr %318, align 8, !alias.scope !34, !noalias !35
  %321 = load double, ptr %319, align 8, !alias.scope !34, !noalias !35
  %322 = insertelement <2 x double> poison, double %320, i64 0
  %323 = insertelement <2 x double> %322, double %321, i64 1
  %324 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %323, <2 x double> %158, <2 x double> %317)
  %325 = getelementptr i8, ptr %194, i64 %45
  %326 = getelementptr i8, ptr %195, i64 %45
  %327 = load double, ptr %325, align 8, !alias.scope !34, !noalias !35
  %328 = load double, ptr %326, align 8, !alias.scope !34, !noalias !35
  %329 = insertelement <2 x double> poison, double %327, i64 0
  %330 = insertelement <2 x double> %329, double %328, i64 1
  %331 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %330, <2 x double> %160, <2 x double> %324)
  %332 = getelementptr i8, ptr %194, i64 %46
  %333 = getelementptr i8, ptr %195, i64 %46
  %334 = load double, ptr %332, align 8, !alias.scope !34, !noalias !35
  %335 = load double, ptr %333, align 8, !alias.scope !34, !noalias !35
  %336 = insertelement <2 x double> poison, double %334, i64 0
  %337 = insertelement <2 x double> %336, double %335, i64 1
  %338 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %337, <2 x double> %162, <2 x double> %331)
  %339 = getelementptr i8, ptr %194, i64 %47
  %340 = getelementptr i8, ptr %195, i64 %47
  %341 = load double, ptr %339, align 8, !alias.scope !34, !noalias !35
  %342 = load double, ptr %340, align 8, !alias.scope !34, !noalias !35
  %343 = insertelement <2 x double> poison, double %341, i64 0
  %344 = insertelement <2 x double> %343, double %342, i64 1
  %345 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %344, <2 x double> %164, <2 x double> %338)
  %346 = getelementptr i8, ptr %194, i64 %48
  %347 = getelementptr i8, ptr %195, i64 %48
  %348 = load double, ptr %346, align 8, !alias.scope !34, !noalias !35
  %349 = load double, ptr %347, align 8, !alias.scope !34, !noalias !35
  %350 = insertelement <2 x double> poison, double %348, i64 0
  %351 = insertelement <2 x double> %350, double %349, i64 1
  %352 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %351, <2 x double> %166, <2 x double> %345)
  %353 = getelementptr i8, ptr %194, i64 %49
  %354 = getelementptr i8, ptr %195, i64 %49
  %355 = load double, ptr %353, align 8, !alias.scope !34, !noalias !35
  %356 = load double, ptr %354, align 8, !alias.scope !34, !noalias !35
  %357 = insertelement <2 x double> poison, double %355, i64 0
  %358 = insertelement <2 x double> %357, double %356, i64 1
  %359 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %358, <2 x double> %168, <2 x double> %352)
  %360 = getelementptr i8, ptr %194, i64 %50
  %361 = getelementptr i8, ptr %195, i64 %50
  %362 = load double, ptr %360, align 8, !alias.scope !34, !noalias !35
  %363 = load double, ptr %361, align 8, !alias.scope !34, !noalias !35
  %364 = insertelement <2 x double> poison, double %362, i64 0
  %365 = insertelement <2 x double> %364, double %363, i64 1
  %366 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %365, <2 x double> %170, <2 x double> %359)
  %367 = getelementptr i8, ptr %194, i64 %51
  %368 = getelementptr i8, ptr %195, i64 %51
  %369 = load double, ptr %367, align 8, !alias.scope !34, !noalias !35
  %370 = load double, ptr %368, align 8, !alias.scope !34, !noalias !35
  %371 = insertelement <2 x double> poison, double %369, i64 0
  %372 = insertelement <2 x double> %371, double %370, i64 1
  %373 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %372, <2 x double> %172, <2 x double> %366)
  %374 = getelementptr i8, ptr %194, i64 %52
  %375 = getelementptr i8, ptr %195, i64 %52
  %376 = load double, ptr %374, align 8, !alias.scope !34, !noalias !35
  %377 = load double, ptr %375, align 8, !alias.scope !34, !noalias !35
  %378 = insertelement <2 x double> poison, double %376, i64 0
  %379 = insertelement <2 x double> %378, double %377, i64 1
  %380 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %379, <2 x double> %174, <2 x double> %373)
  %381 = getelementptr i8, ptr %194, i64 %53
  %382 = getelementptr i8, ptr %195, i64 %53
  %383 = load double, ptr %381, align 8, !alias.scope !34, !noalias !35
  %384 = load double, ptr %382, align 8, !alias.scope !34, !noalias !35
  %385 = insertelement <2 x double> poison, double %383, i64 0
  %386 = insertelement <2 x double> %385, double %384, i64 1
  %387 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %386, <2 x double> %176, <2 x double> %380)
  %388 = getelementptr i8, ptr %194, i64 %54
  %389 = getelementptr i8, ptr %195, i64 %54
  %390 = load double, ptr %388, align 8, !alias.scope !34, !noalias !35
  %391 = load double, ptr %389, align 8, !alias.scope !34, !noalias !35
  %392 = insertelement <2 x double> poison, double %390, i64 0
  %393 = insertelement <2 x double> %392, double %391, i64 1
  %394 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %393, <2 x double> %178, <2 x double> %387)
  %395 = getelementptr i8, ptr %194, i64 %55
  %396 = getelementptr i8, ptr %195, i64 %55
  %397 = load double, ptr %395, align 8, !alias.scope !34, !noalias !35
  %398 = load double, ptr %396, align 8, !alias.scope !34, !noalias !35
  %399 = insertelement <2 x double> poison, double %397, i64 0
  %400 = insertelement <2 x double> %399, double %398, i64 1
  %401 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %400, <2 x double> %180, <2 x double> %394)
  %402 = getelementptr i8, ptr %194, i64 %56
  %403 = getelementptr i8, ptr %195, i64 %56
  %404 = load double, ptr %402, align 8, !alias.scope !34, !noalias !35
  %405 = load double, ptr %403, align 8, !alias.scope !34, !noalias !35
  %406 = insertelement <2 x double> poison, double %404, i64 0
  %407 = insertelement <2 x double> %406, double %405, i64 1
  %408 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %407, <2 x double> %182, <2 x double> %401)
  %409 = getelementptr i8, ptr %194, i64 %57
  %410 = getelementptr i8, ptr %195, i64 %57
  %411 = load double, ptr %409, align 8, !alias.scope !34, !noalias !35
  %412 = load double, ptr %410, align 8, !alias.scope !34, !noalias !35
  %413 = insertelement <2 x double> poison, double %411, i64 0
  %414 = insertelement <2 x double> %413, double %412, i64 1
  %415 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %414, <2 x double> %184, <2 x double> %408)
  %416 = getelementptr i8, ptr %194, i64 %58
  %417 = getelementptr i8, ptr %195, i64 %58
  %418 = load double, ptr %416, align 8, !alias.scope !34, !noalias !35
  %419 = load double, ptr %417, align 8, !alias.scope !34, !noalias !35
  %420 = insertelement <2 x double> poison, double %418, i64 0
  %421 = insertelement <2 x double> %420, double %419, i64 1
  %422 = call <2 x double> @llvm.fmuladd.v2f64(<2 x double> %421, <2 x double> %186, <2 x double> %415)
  store <2 x double> %422, ptr %197, align 8, !alias.scope !22, !noalias !23
  %423 = add nuw i64 %188, 2
  %424 = icmp eq i64 %423, 32
  br i1 %424, label %425, label %187, !llvm.loop !36

425:                                              ; preds = %187
  %426 = add nuw nsw i64 %26, 1
  %427 = icmp eq i64 %426, 625
  br i1 %427, label %15, label %25
}

; Function Attrs: nocallback nofree nounwind willreturn memory(argmem: write)
declare void @llvm.memset.p0.i64(ptr nocapture writeonly, i8, i64, i1 immarg) #7

; Function Attrs: nocallback nofree nosync nounwind speculatable willreturn memory(none)
declare <2 x double> @llvm.fmuladd.v2f64(<2 x double>, <2 x double>, <2 x double>) #8

attributes #0 = { nounwind uwtable "min-legal-vector-width"="0" "no-trapping-math"="true" "polly-optimized" "stack-protector-buffer-size"="8" "target-cpu"="x86-64" "target-features"="+cmov,+cx8,+fxsr,+mmx,+sse,+sse2,+x87" "tune-cpu"="generic" }
attributes #1 = { "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="x86-64" "target-features"="+cmov,+cx8,+fxsr,+mmx,+sse,+sse2,+x87" "tune-cpu"="generic" }
attributes #2 = { mustprogress nounwind willreturn allockind("free") memory(argmem: readwrite, inaccessiblemem: readwrite) "alloc-family"="malloc" "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="x86-64" "target-features"="+cmov,+cx8,+fxsr,+mmx,+sse,+sse2,+x87" "tune-cpu"="generic" }
attributes #3 = { mustprogress nocallback nofree nosync nounwind speculatable willreturn memory(none) }
attributes #4 = { nofree nounwind "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="x86-64" "target-features"="+cmov,+cx8,+fxsr,+mmx,+sse,+sse2,+x87" "tune-cpu"="generic" }
attributes #5 = { nofree nounwind }
attributes #6 = { "polly.skip.fn" }
attributes #7 = { nocallback nofree nounwind willreturn memory(argmem: write) }
attributes #8 = { nocallback nofree nosync nounwind speculatable willreturn memory(none) }
attributes #9 = { nounwind }
attributes #10 = { cold }

!llvm.module.flags = !{!0, !1, !2, !3, !4}
!llvm.ident = !{!5}

!0 = !{i32 1, !"wchar_size", i32 4}
!1 = !{i32 7, !"openmp", i32 51}
!2 = !{i32 8, !"PIC Level", i32 2}
!3 = !{i32 7, !"PIE Level", i32 2}
!4 = !{i32 7, !"uwtable", i32 2}
!5 = !{!"clang version 18.1.3 (https://github.com/llvm/llvm-project.git c13b7485b87909fcf739f62cfa382b55407433c0)"}
!6 = !{!7, !7, i64 0}
!7 = !{!"double", !8, i64 0}
!8 = !{!"omnipotent char", !9, i64 0}
!9 = !{!"Simple C/C++ TBAA"}
!10 = distinct !{!10, !11, !12}
!11 = !{!"llvm.loop.isvectorized", i32 1}
!12 = !{!"llvm.loop.unroll.runtime.disable"}
!13 = distinct !{!13, !11}
!14 = !{!15}
!15 = distinct !{!15, !16, !"polly.alias.scope.MemRef1"}
!16 = distinct !{!16, !"polly.alias.scope.domain"}
!17 = !{!18, !19, !20, !21}
!18 = distinct !{!18, !16, !"polly.alias.scope.MemRef0"}
!19 = distinct !{!19, !16, !"polly.alias.scope.MemRef2"}
!20 = distinct !{!20, !16, !"polly.alias.scope.MemRef3"}
!21 = distinct !{!21, !16, !"polly.alias.scope.MemRef4"}
!22 = !{!18}
!23 = !{!15, !19, !20, !21}
!24 = distinct !{}
!25 = distinct !{!25, !26, !11, !12}
!26 = !{!"llvm.loop.parallel_accesses", !24}
!27 = !{!28, !28, i64 0}
!28 = !{!"any pointer", !8, i64 0}
!29 = !{!20}
!30 = !{!18, !15, !19, !21}
!31 = !{!21}
!32 = !{!18, !15, !19, !20}
!33 = distinct !{!33, !11, !12}
!34 = !{!19}
!35 = !{!18, !15, !20, !21}
!36 = distinct !{!36, !11, !12}
