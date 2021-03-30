; ModuleID = 'enzymeSphericalHarmonics.cpp'
source_filename = "enzymeSphericalHarmonics.cpp"
target datalayout = "e-m:o-p270:32:32-p271:32:32-p272:64:64-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-apple-macosx10.15.0"

@.str = private unnamed_addr constant [12 x i8] c"Ylmr[0]=%f\0A\00", align 1
@enzyme_const = external local_unnamed_addr global i32, align 4
@.str.1 = private unnamed_addr constant [14 x i8] c"d_Ylmr[0]=%f\0A\00", align 1
@.memset_pattern = private unnamed_addr constant [2 x double] [double 0x3FD20DD750429B6D, double 0x3FD20DD750429B6D], align 16
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z21cpuSphericalHarmonicsIdEvPT_S1_S1_S1_S1_S1_S1_S1_S0_ii(double* %0, double* %1, double* %2, double* %3, double* %4, double* %5, double* %6, double* %7, double %8, i32 %9, i32 %10) #0 {
  %12 = icmp sgt i32 %10, 0
  br i1 %12, label %13, label %27

13:                                               ; preds = %11
  %14 = fmul fast double %8, 4.000000e+00
  %15 = fdiv fast double 2.500000e-01, %8
  %16 = tail call fast double @llvm.sqrt.f64(double %15)
  %17 = getelementptr inbounds double, double* %4, i64 1
  %18 = getelementptr inbounds double, double* %6, i64 1
  %19 = getelementptr inbounds double, double* %6, i64 2
  %20 = getelementptr inbounds double, double* %7, i64 1
  %21 = shl nuw i32 %10, 1
  %22 = icmp slt i32 %9, 2
  %23 = zext i32 %10 to i64
  %24 = add i32 %9, 1
  %25 = sext i32 %21 to i64
  %26 = zext i32 %24 to i64
  br label %28

27:                                               ; preds = %161, %11
  ret void

28:                                               ; preds = %161, %13
  %29 = phi i64 [ 0, %13 ], [ %162, %161 ]
  %30 = getelementptr inbounds double, double* %0, i64 %29
  store double %16, double* %30, align 8, !tbaa !4
  %31 = getelementptr inbounds double, double* %1, i64 %29
  store double 0.000000e+00, double* %31, align 8, !tbaa !4
  %32 = getelementptr inbounds double, double* %2, i64 %29
  %33 = load double, double* %32, align 8, !tbaa !4
  %34 = tail call fast double @llvm.cos.f64(double %33)
  %35 = tail call fast double @llvm.sin.f64(double %33)
  %36 = fneg fast double %35
  store double %34, double* %4, align 8, !tbaa !4
  store double %36, double* %17, align 8, !tbaa !4
  %37 = load double, double* %18, align 8, !tbaa !4
  %38 = fmul fast double %37, 3.000000e+00
  %39 = fmul fast double %37, %14
  %40 = fdiv fast double %38, %39
  %41 = tail call fast double @llvm.sqrt.f64(double %40)
  store double %41, double* %7, align 8, !tbaa !4
  %42 = load double, double* %4, align 8, !tbaa !4
  %43 = fmul fast double %41, %42
  %44 = add nuw nsw i64 %29, %23
  %45 = getelementptr inbounds double, double* %0, i64 %44
  store double %43, double* %45, align 8, !tbaa !4
  %46 = getelementptr inbounds double, double* %1, i64 %44
  store double 0.000000e+00, double* %46, align 8, !tbaa !4
  %47 = load double, double* %6, align 8, !tbaa !4
  %48 = fmul fast double %47, 3.000000e+00
  %49 = load double, double* %19, align 8, !tbaa !4
  %50 = fmul fast double %49, %14
  %51 = fdiv fast double %48, %50
  %52 = tail call fast double @llvm.sqrt.f64(double %51)
  store double %52, double* %20, align 8, !tbaa !4
  %53 = getelementptr inbounds double, double* %3, i64 %29
  %54 = load double, double* %53, align 8, !tbaa !4
  %55 = tail call fast double @llvm.cos.f64(double %54)
  %56 = fmul fast double %52, %55
  %57 = load double, double* %17, align 8, !tbaa !4
  %58 = fmul fast double %56, %57
  %59 = add nsw i64 %29, %25
  %60 = getelementptr inbounds double, double* %0, i64 %59
  store double %58, double* %60, align 8, !tbaa !4
  %61 = load double, double* %20, align 8, !tbaa !4
  %62 = load double, double* %53, align 8, !tbaa !4
  %63 = tail call fast double @llvm.sin.f64(double %62)
  %64 = fmul fast double %63, %61
  %65 = load double, double* %17, align 8, !tbaa !4
  %66 = fmul fast double %64, %65
  %67 = getelementptr inbounds double, double* %1, i64 %59
  store double %66, double* %67, align 8, !tbaa !4
  br i1 %22, label %161, label %68

68:                                               ; preds = %28
  %69 = insertelement <2 x double> undef, double %34, i32 0
  br label %70

70:                                               ; preds = %68, %157
  %71 = phi i64 [ %115, %157 ], [ 2, %68 ]
  %72 = phi i64 [ %158, %157 ], [ 1, %68 ]
  %73 = phi i64 [ %159, %157 ], [ 3, %68 ]
  %74 = phi i32 [ %116, %157 ], [ 2, %68 ]
  %75 = add nsw i64 %71, -1
  %76 = getelementptr inbounds double, double* %4, i64 %75
  %77 = load double, double* %76, align 8, !tbaa !4
  %78 = getelementptr inbounds double, double* %5, i64 %75
  store double %77, double* %78, align 8, !tbaa !4
  %79 = trunc i64 %75 to i32
  %80 = shl i32 %79, 1
  %81 = or i32 %80, 1
  %82 = sitofp i32 %81 to double
  %83 = insertelement <2 x double> %69, double %77, i32 1
  %84 = insertelement <2 x double> undef, double %82, i32 0
  %85 = shufflevector <2 x double> %84, <2 x double> undef, <2 x i32> zeroinitializer
  %86 = fmul fast <2 x double> %83, %85
  %87 = insertelement <2 x double> undef, double %77, i32 0
  %88 = insertelement <2 x double> %87, double %36, i32 1
  %89 = fmul fast <2 x double> %86, %88
  %90 = bitcast double* %76 to <2 x double>*
  store <2 x double> %89, <2 x double>* %90, align 8, !tbaa !4
  %91 = extractelement <2 x double> %86, i32 0
  br label %92

92:                                               ; preds = %92, %70
  %93 = phi i64 [ 0, %70 ], [ %108, %92 ]
  %94 = getelementptr inbounds double, double* %4, i64 %93
  %95 = load double, double* %94, align 8, !tbaa !4
  %96 = fmul fast double %91, %95
  %97 = add nuw nsw i64 %75, %93
  %98 = trunc i64 %97 to i32
  %99 = sitofp i32 %98 to double
  %100 = getelementptr inbounds double, double* %5, i64 %93
  %101 = load double, double* %100, align 8, !tbaa !4
  %102 = fmul fast double %101, %99
  %103 = fsub fast double %96, %102
  %104 = sub nsw i64 %71, %93
  %105 = trunc i64 %104 to i32
  %106 = sitofp i32 %105 to double
  %107 = fdiv fast double %103, %106
  store double %107, double* %94, align 8, !tbaa !4
  store double %95, double* %100, align 8, !tbaa !4
  %108 = add nuw nsw i64 %93, 1
  %109 = icmp eq i64 %108, %72
  br i1 %109, label %110, label %92

110:                                              ; preds = %92
  %111 = trunc i64 %71 to i32
  %112 = shl i32 %111, 1
  %113 = or i32 %112, 1
  %114 = sitofp i32 %113 to double
  %115 = add nuw nsw i64 %71, 1
  %116 = add nuw nsw i32 %74, 1
  %117 = mul nsw i32 %116, %111
  %118 = lshr i32 %117, 1
  %119 = zext i32 %118 to i64
  br label %120

120:                                              ; preds = %120, %110
  %121 = phi i64 [ %155, %120 ], [ 0, %110 ]
  %122 = sub nsw i64 %71, %121
  %123 = getelementptr inbounds double, double* %6, i64 %122
  %124 = load double, double* %123, align 8, !tbaa !4
  %125 = fmul fast double %124, %114
  %126 = add nuw nsw i64 %121, %71
  %127 = and i64 %126, 4294967295
  %128 = getelementptr inbounds double, double* %6, i64 %127
  %129 = load double, double* %128, align 8, !tbaa !4
  %130 = fmul fast double %129, %14
  %131 = fdiv fast double %125, %130
  %132 = tail call fast double @llvm.sqrt.f64(double %131)
  %133 = getelementptr inbounds double, double* %7, i64 %121
  store double %132, double* %133, align 8, !tbaa !4
  %134 = trunc i64 %121 to i32
  %135 = sitofp i32 %134 to double
  %136 = load double, double* %53, align 8, !tbaa !4
  %137 = fmul fast double %136, %135
  %138 = tail call fast double @llvm.cos.f64(double %137)
  %139 = getelementptr inbounds double, double* %4, i64 %121
  %140 = load double, double* %139, align 8, !tbaa !4
  %141 = fmul fast double %132, %140
  %142 = fmul fast double %141, %138
  %143 = add nuw nsw i64 %121, %119
  %144 = mul nsw i64 %143, %23
  %145 = add nuw nsw i64 %144, %29
  %146 = getelementptr inbounds double, double* %0, i64 %145
  store double %142, double* %146, align 8, !tbaa !4
  %147 = load double, double* %133, align 8, !tbaa !4
  %148 = load double, double* %53, align 8, !tbaa !4
  %149 = fmul fast double %148, %135
  %150 = tail call fast double @llvm.sin.f64(double %149)
  %151 = load double, double* %139, align 8, !tbaa !4
  %152 = fmul fast double %151, %147
  %153 = fmul fast double %152, %150
  %154 = getelementptr inbounds double, double* %1, i64 %145
  store double %153, double* %154, align 8, !tbaa !4
  %155 = add nuw nsw i64 %121, 1
  %156 = icmp eq i64 %155, %73
  br i1 %156, label %157, label %120

157:                                              ; preds = %120
  %158 = add nuw nsw i64 %72, 1
  %159 = add nuw nsw i64 %73, 1
  %160 = icmp eq i64 %115, %26
  br i1 %160, label %161, label %70

161:                                              ; preds = %157, %28
  %162 = add nuw nsw i64 %29, 1
  %163 = icmp eq i64 %162, %23
  br i1 %163, label %27, label %28
}
; Function Attrs: argmemonly nounwind willreturn
declare void @llvm.lifetime.start.p0i8(i64 immarg, i8* nocapture) #1
; Function Attrs: nounwind readnone speculatable willreturn
declare double @llvm.sqrt.f64(double) #2
; Function Attrs: nounwind readnone speculatable willreturn
declare double @llvm.cos.f64(double) #2
; Function Attrs: nounwind readnone speculatable willreturn
declare double @llvm.sin.f64(double) #2
; Function Attrs: argmemonly nounwind willreturn
declare void @llvm.lifetime.end.p0i8(i64 immarg, i8* nocapture) #1
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z21cpuSphericalHarmonicsIfEvPT_S1_S1_S1_S1_S1_S1_S1_S0_ii(float* %0, float* %1, float* %2, float* %3, float* %4, float* %5, float* %6, float* %7, float %8, i32 %9, i32 %10) local_unnamed_addr #0 {
  %12 = icmp sgt i32 %10, 0
  br i1 %12, label %13, label %27

13:                                               ; preds = %11
  %14 = fmul fast float %8, 4.000000e+00
  %15 = fdiv fast float 2.500000e-01, %8
  %16 = tail call fast float @llvm.sqrt.f32(float %15) #8
  %17 = getelementptr inbounds float, float* %4, i64 1
  %18 = getelementptr inbounds float, float* %6, i64 1
  %19 = getelementptr inbounds float, float* %6, i64 2
  %20 = getelementptr inbounds float, float* %7, i64 1
  %21 = shl nuw i32 %10, 1
  %22 = icmp slt i32 %9, 2
  %23 = zext i32 %10 to i64
  %24 = add i32 %9, 1
  %25 = sext i32 %21 to i64
  %26 = zext i32 %24 to i64
  br label %28

27:                                               ; preds = %202, %11
  ret void

28:                                               ; preds = %202, %13
  %29 = phi i64 [ 0, %13 ], [ %203, %202 ]
  %30 = getelementptr inbounds float, float* %0, i64 %29
  store float %16, float* %30, align 4, !tbaa !8
  %31 = getelementptr inbounds float, float* %1, i64 %29
  store float 0.000000e+00, float* %31, align 4, !tbaa !8
  %32 = getelementptr inbounds float, float* %2, i64 %29
  %33 = load float, float* %32, align 4, !tbaa !8
  %34 = tail call fast float @llvm.cos.f32(float %33) #8
  %35 = tail call fast float @llvm.sin.f32(float %33) #8
  %36 = fneg fast float %35
  store float %34, float* %4, align 4, !tbaa !8
  store float %36, float* %17, align 4, !tbaa !8
  %37 = load float, float* %18, align 4, !tbaa !8
  %38 = fmul fast float %37, 3.000000e+00
  %39 = fmul fast float %37, %14
  %40 = fdiv fast float %38, %39
  %41 = tail call fast float @llvm.sqrt.f32(float %40) #8
  store float %41, float* %7, align 4, !tbaa !8
  %42 = load float, float* %4, align 4, !tbaa !8
  %43 = fmul fast float %41, %42
  %44 = add nuw nsw i64 %29, %23
  %45 = getelementptr inbounds float, float* %0, i64 %44
  store float %43, float* %45, align 4, !tbaa !8
  %46 = getelementptr inbounds float, float* %1, i64 %44
  store float 0.000000e+00, float* %46, align 4, !tbaa !8
  %47 = load float, float* %6, align 4, !tbaa !8
  %48 = fmul fast float %47, 3.000000e+00
  %49 = load float, float* %19, align 4, !tbaa !8
  %50 = fmul fast float %49, %14
  %51 = fdiv fast float %48, %50
  %52 = tail call fast float @llvm.sqrt.f32(float %51) #8
  store float %52, float* %20, align 4, !tbaa !8
  %53 = getelementptr inbounds float, float* %3, i64 %29
  %54 = load float, float* %53, align 4, !tbaa !8
  %55 = tail call fast float @llvm.cos.f32(float %54) #8
  %56 = fmul fast float %52, %55
  %57 = load float, float* %17, align 4, !tbaa !8
  %58 = fmul fast float %56, %57
  %59 = add nsw i64 %29, %25
  %60 = getelementptr inbounds float, float* %0, i64 %59
  store float %58, float* %60, align 4, !tbaa !8
  %61 = load float, float* %20, align 4, !tbaa !8
  %62 = load float, float* %53, align 4, !tbaa !8
  %63 = tail call fast float @llvm.sin.f32(float %62) #8
  %64 = fmul fast float %63, %61
  %65 = load float, float* %17, align 4, !tbaa !8
  %66 = fmul fast float %64, %65
  %67 = getelementptr inbounds float, float* %1, i64 %59
  store float %66, float* %67, align 4, !tbaa !8
  br i1 %22, label %202, label %68

68:                                               ; preds = %28, %197
  %69 = phi i64 [ %201, %197 ], [ 0, %28 ]
  %70 = phi i64 [ %155, %197 ], [ 2, %28 ]
  %71 = phi i64 [ %198, %197 ], [ 1, %28 ]
  %72 = phi i64 [ %199, %197 ], [ 3, %28 ]
  %73 = phi i32 [ %156, %197 ], [ 2, %28 ]
  %74 = add i64 %69, 1
  %75 = add i64 %69, 1
  %76 = add nsw i64 %70, -1
  %77 = getelementptr inbounds float, float* %4, i64 %76
  %78 = load float, float* %77, align 4, !tbaa !8
  %79 = getelementptr inbounds float, float* %5, i64 %76
  store float %78, float* %79, align 4, !tbaa !8
  %80 = trunc i64 %76 to i32
  %81 = shl i32 %80, 1
  %82 = or i32 %81, 1
  %83 = sitofp i32 %82 to float
  %84 = fmul fast float %34, %83
  %85 = fmul fast float %84, %78
  store float %85, float* %77, align 4, !tbaa !8
  %86 = fmul fast float %78, %83
  %87 = fmul fast float %86, %36
  %88 = getelementptr inbounds float, float* %4, i64 %70
  store float %87, float* %88, align 4, !tbaa !8
  %89 = icmp ult i64 %75, 4
  br i1 %89, label %130, label %90

90:                                               ; preds = %68
  %91 = getelementptr float, float* %5, i64 %74
  %92 = getelementptr float, float* %4, i64 %74
  %93 = icmp ugt float* %91, %4
  %94 = icmp ugt float* %92, %5
  %95 = and i1 %93, %94
  br i1 %95, label %130, label %96

96:                                               ; preds = %90
  %97 = and i64 %75, -4
  %98 = insertelement <4 x float> undef, float %84, i32 0
  %99 = shufflevector <4 x float> %98, <4 x float> undef, <4 x i32> zeroinitializer
  %100 = insertelement <4 x i64> undef, i64 %76, i32 0
  %101 = shufflevector <4 x i64> %100, <4 x i64> undef, <4 x i32> zeroinitializer
  %102 = insertelement <4 x i64> undef, i64 %70, i32 0
  %103 = shufflevector <4 x i64> %102, <4 x i64> undef, <4 x i32> zeroinitializer
  br label %104

104:                                              ; preds = %104, %96
  %105 = phi i64 [ 0, %96 ], [ %125, %104 ]
  %106 = phi <4 x i64> [ <i64 0, i64 1, i64 2, i64 3>, %96 ], [ %126, %104 ]
  %107 = getelementptr inbounds float, float* %4, i64 %105
  %108 = bitcast float* %107 to <4 x float>*
  %109 = load <4 x float>, <4 x float>* %108, align 4, !tbaa !8, !alias.scope !10, !noalias !13
  %110 = fmul fast <4 x float> %99, %109
  %111 = add nuw nsw <4 x i64> %101, %106
  %112 = trunc <4 x i64> %111 to <4 x i32>
  %113 = sitofp <4 x i32> %112 to <4 x float>
  %114 = getelementptr inbounds float, float* %5, i64 %105
  %115 = bitcast float* %114 to <4 x float>*
  %116 = load <4 x float>, <4 x float>* %115, align 4, !tbaa !8, !alias.scope !13
  %117 = fmul fast <4 x float> %116, %113
  %118 = fsub fast <4 x float> %110, %117
  %119 = sub nsw <4 x i64> %103, %106
  %120 = trunc <4 x i64> %119 to <4 x i32>
  %121 = sitofp <4 x i32> %120 to <4 x float>
  %122 = fdiv fast <4 x float> %118, %121
  %123 = bitcast float* %107 to <4 x float>*
  store <4 x float> %122, <4 x float>* %123, align 4, !tbaa !8, !alias.scope !10, !noalias !13
  %124 = bitcast float* %114 to <4 x float>*
  store <4 x float> %109, <4 x float>* %124, align 4, !tbaa !8, !alias.scope !13
  %125 = add i64 %105, 4
  %126 = add <4 x i64> %106, <i64 4, i64 4, i64 4, i64 4>
  %127 = icmp eq i64 %125, %97
  br i1 %127, label %128, label %104, !llvm.loop !15

128:                                              ; preds = %104
  %129 = icmp eq i64 %75, %97
  br i1 %129, label %150, label %130

130:                                              ; preds = %128, %90, %68
  %131 = phi i64 [ 0, %90 ], [ 0, %68 ], [ %97, %128 ]
  br label %132

132:                                              ; preds = %130, %132
  %133 = phi i64 [ %148, %132 ], [ %131, %130 ]
  %134 = getelementptr inbounds float, float* %4, i64 %133
  %135 = load float, float* %134, align 4, !tbaa !8
  %136 = fmul fast float %84, %135
  %137 = add nuw nsw i64 %76, %133
  %138 = trunc i64 %137 to i32
  %139 = sitofp i32 %138 to float
  %140 = getelementptr inbounds float, float* %5, i64 %133
  %141 = load float, float* %140, align 4, !tbaa !8
  %142 = fmul fast float %141, %139
  %143 = fsub fast float %136, %142
  %144 = sub nsw i64 %70, %133
  %145 = trunc i64 %144 to i32
  %146 = sitofp i32 %145 to float
  %147 = fdiv fast float %143, %146
  store float %147, float* %134, align 4, !tbaa !8
  store float %135, float* %140, align 4, !tbaa !8
  %148 = add nuw nsw i64 %133, 1
  %149 = icmp eq i64 %148, %71
  br i1 %149, label %150, label %132, !llvm.loop !17

150:                                              ; preds = %132, %128
  %151 = trunc i64 %70 to i32
  %152 = shl i32 %151, 1
  %153 = or i32 %152, 1
  %154 = sitofp i32 %153 to float
  %155 = add nuw nsw i64 %70, 1
  %156 = add nuw nsw i32 %73, 1
  %157 = mul nsw i32 %156, %151
  %158 = lshr i32 %157, 1
  %159 = zext i32 %158 to i64
  br label %160

160:                                              ; preds = %160, %150
  %161 = phi i64 [ %195, %160 ], [ 0, %150 ]
  %162 = sub nsw i64 %70, %161
  %163 = getelementptr inbounds float, float* %6, i64 %162
  %164 = load float, float* %163, align 4, !tbaa !8
  %165 = fmul fast float %164, %154
  %166 = add nuw nsw i64 %161, %70
  %167 = and i64 %166, 4294967295
  %168 = getelementptr inbounds float, float* %6, i64 %167
  %169 = load float, float* %168, align 4, !tbaa !8
  %170 = fmul fast float %169, %14
  %171 = fdiv fast float %165, %170
  %172 = tail call fast float @llvm.sqrt.f32(float %171) #8
  %173 = getelementptr inbounds float, float* %7, i64 %161
  store float %172, float* %173, align 4, !tbaa !8
  %174 = trunc i64 %161 to i32
  %175 = sitofp i32 %174 to float
  %176 = load float, float* %53, align 4, !tbaa !8
  %177 = fmul fast float %176, %175
  %178 = tail call fast float @llvm.cos.f32(float %177) #8
  %179 = getelementptr inbounds float, float* %4, i64 %161
  %180 = load float, float* %179, align 4, !tbaa !8
  %181 = fmul fast float %172, %180
  %182 = fmul fast float %181, %178
  %183 = add nuw nsw i64 %161, %159
  %184 = mul nsw i64 %183, %23
  %185 = add nuw nsw i64 %184, %29
  %186 = getelementptr inbounds float, float* %0, i64 %185
  store float %182, float* %186, align 4, !tbaa !8
  %187 = load float, float* %173, align 4, !tbaa !8
  %188 = load float, float* %53, align 4, !tbaa !8
  %189 = fmul fast float %188, %175
  %190 = tail call fast float @llvm.sin.f32(float %189) #8
  %191 = load float, float* %179, align 4, !tbaa !8
  %192 = fmul fast float %191, %187
  %193 = fmul fast float %192, %190
  %194 = getelementptr inbounds float, float* %1, i64 %185
  store float %193, float* %194, align 4, !tbaa !8
  %195 = add nuw nsw i64 %161, 1
  %196 = icmp eq i64 %195, %72
  br i1 %196, label %197, label %160

197:                                              ; preds = %160
  %198 = add nuw nsw i64 %71, 1
  %199 = add nuw nsw i64 %72, 1
  %200 = icmp eq i64 %155, %26
  %201 = add i64 %69, 1
  br i1 %200, label %202, label %68

202:                                              ; preds = %197, %28
  %203 = add nuw nsw i64 %29, 1
  %204 = icmp eq i64 %203, %23
  br i1 %204, label %27, label %28
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z26cpuSphericalHarmonicsDerivIdEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S0_ii(double* %0, double* %1, double* %2, double* %3, double* %4, double* %5, double* %6, double* %7, double* %8, double* %9, double* %10, double* %11, double* %12, double* %13, double %14, i32 %15, i32 %16) local_unnamed_addr #0 {
  %18 = icmp sgt i32 %16, 0
  br i1 %18, label %19, label %34

19:                                               ; preds = %17
  %20 = fmul fast double %14, 4.000000e+00
  %21 = fdiv fast double 2.500000e-01, %14
  %22 = tail call fast double @llvm.sqrt.f64(double %21)
  %23 = getelementptr inbounds double, double* %8, i64 1
  %24 = getelementptr inbounds double, double* %12, i64 1
  %25 = getelementptr inbounds double, double* %12, i64 2
  %26 = getelementptr inbounds double, double* %13, i64 1
  %27 = shl nuw i32 %16, 1
  %28 = icmp slt i32 %15, 2
  %29 = add i32 %15, 1
  %30 = sext i32 %27 to i64
  %31 = zext i32 %16 to i64
  %32 = zext i32 %29 to i64
  %33 = bitcast double* %10 to <2 x double>*
  br label %35

34:                                               ; preds = %357, %17
  ret void

35:                                               ; preds = %357, %19
  %36 = phi i64 [ 0, %19 ], [ %358, %357 ]
  %37 = getelementptr inbounds double, double* %0, i64 %36
  store double %22, double* %37, align 8, !tbaa !4
  %38 = getelementptr inbounds double, double* %1, i64 %36
  store double 0.000000e+00, double* %38, align 8, !tbaa !4
  %39 = getelementptr inbounds double, double* %2, i64 %36
  store double 0.000000e+00, double* %39, align 8, !tbaa !4
  %40 = getelementptr inbounds double, double* %4, i64 %36
  store double 0.000000e+00, double* %40, align 8, !tbaa !4
  %41 = getelementptr inbounds double, double* %3, i64 %36
  store double 0.000000e+00, double* %41, align 8, !tbaa !4
  %42 = getelementptr inbounds double, double* %5, i64 %36
  store double 0.000000e+00, double* %42, align 8, !tbaa !4
  %43 = getelementptr inbounds double, double* %6, i64 %36
  %44 = load double, double* %43, align 8, !tbaa !4
  %45 = tail call fast double @llvm.cos.f64(double %44)
  %46 = tail call fast double @llvm.sin.f64(double %44)
  %47 = insertelement <2 x double> undef, double %46, i32 0
  %48 = insertelement <2 x double> %47, double %45, i32 1
  %49 = fneg fast <2 x double> %48
  store double %45, double* %8, align 8, !tbaa !4
  %50 = extractelement <2 x double> %49, i32 0
  store double %50, double* %23, align 8, !tbaa !4
  store <2 x double> %49, <2 x double>* %33, align 8, !tbaa !4
  %51 = load double, double* %24, align 8, !tbaa !4
  %52 = fmul fast double %51, 3.000000e+00
  %53 = fmul fast double %51, %20
  %54 = fdiv fast double %52, %53
  %55 = tail call fast double @llvm.sqrt.f64(double %54)
  store double %55, double* %13, align 8, !tbaa !4
  %56 = load double, double* %8, align 8, !tbaa !4
  %57 = fmul fast double %55, %56
  %58 = add nuw nsw i64 %36, %31
  %59 = getelementptr inbounds double, double* %0, i64 %58
  store double %57, double* %59, align 8, !tbaa !4
  %60 = getelementptr inbounds double, double* %1, i64 %58
  store double 0.000000e+00, double* %60, align 8, !tbaa !4
  %61 = load double, double* %13, align 8, !tbaa !4
  %62 = fmul fast double %61, %50
  %63 = getelementptr inbounds double, double* %2, i64 %58
  store double %62, double* %63, align 8, !tbaa !4
  %64 = getelementptr inbounds double, double* %3, i64 %58
  store double 0.000000e+00, double* %64, align 8, !tbaa !4
  %65 = getelementptr inbounds double, double* %4, i64 %58
  store double 0.000000e+00, double* %65, align 8, !tbaa !4
  %66 = getelementptr inbounds double, double* %5, i64 %58
  store double 0.000000e+00, double* %66, align 8, !tbaa !4
  %67 = load double, double* %12, align 8, !tbaa !4
  %68 = fmul fast double %67, 3.000000e+00
  %69 = load double, double* %25, align 8, !tbaa !4
  %70 = fmul fast double %69, %20
  %71 = fdiv fast double %68, %70
  %72 = tail call fast double @llvm.sqrt.f64(double %71)
  store double %72, double* %26, align 8, !tbaa !4
  %73 = getelementptr inbounds double, double* %7, i64 %36
  %74 = load double, double* %73, align 8, !tbaa !4
  %75 = tail call fast double @llvm.cos.f64(double %74)
  %76 = fmul fast double %72, %75
  %77 = load double, double* %23, align 8, !tbaa !4
  %78 = fmul fast double %76, %77
  %79 = add nsw i64 %36, %30
  %80 = getelementptr inbounds double, double* %0, i64 %79
  store double %78, double* %80, align 8, !tbaa !4
  %81 = load double, double* %26, align 8, !tbaa !4
  %82 = load double, double* %73, align 8, !tbaa !4
  %83 = tail call fast double @llvm.sin.f64(double %82)
  %84 = fmul fast double %83, %81
  %85 = load double, double* %23, align 8, !tbaa !4
  %86 = fmul fast double %84, %85
  %87 = getelementptr inbounds double, double* %1, i64 %79
  store double %86, double* %87, align 8, !tbaa !4
  %88 = load double, double* %26, align 8, !tbaa !4
  %89 = load double, double* %73, align 8, !tbaa !4
  %90 = tail call fast double @llvm.cos.f64(double %89)
  %91 = extractelement <2 x double> %49, i32 1
  %92 = fmul fast double %88, %91
  %93 = fmul fast double %92, %90
  %94 = getelementptr inbounds double, double* %2, i64 %79
  store double %93, double* %94, align 8, !tbaa !4
  %95 = load double, double* %26, align 8, !tbaa !4
  %96 = load double, double* %73, align 8, !tbaa !4
  %97 = tail call fast double @llvm.sin.f64(double %96)
  %98 = fmul fast double %95, %91
  %99 = fmul fast double %98, %97
  %100 = getelementptr inbounds double, double* %3, i64 %79
  store double %99, double* %100, align 8, !tbaa !4
  %101 = load double, double* %26, align 8, !tbaa !4
  %102 = load double, double* %73, align 8, !tbaa !4
  %103 = tail call fast double @llvm.sin.f64(double %102)
  %104 = load double, double* %23, align 8, !tbaa !4
  %105 = fsub fast double -0.000000e+00, %101
  %106 = fmul fast double %104, %105
  %107 = fmul fast double %106, %103
  %108 = getelementptr inbounds double, double* %4, i64 %79
  store double %107, double* %108, align 8, !tbaa !4
  %109 = load double, double* %26, align 8, !tbaa !4
  %110 = load double, double* %73, align 8, !tbaa !4
  %111 = tail call fast double @llvm.cos.f64(double %110)
  %112 = load double, double* %23, align 8, !tbaa !4
  %113 = fmul fast double %112, %109
  %114 = fmul fast double %113, %111
  %115 = getelementptr inbounds double, double* %5, i64 %79
  store double %114, double* %115, align 8, !tbaa !4
  br i1 %28, label %357, label %116

116:                                              ; preds = %35
  %117 = trunc i64 %36 to i32
  %118 = insertelement <2 x double> undef, double %45, i32 0
  %119 = insertelement <2 x double> undef, double %45, i32 0
  %120 = shufflevector <2 x double> %119, <2 x double> undef, <2 x i32> zeroinitializer
  %121 = shufflevector <2 x double> %47, <2 x double> undef, <2 x i32> zeroinitializer
  br label %127

122:                                              ; preds = %282
  %123 = add nuw nsw i64 %130, 1
  %124 = add nuw nsw i64 %131, 1
  %125 = icmp eq i64 %274, %32
  %126 = add i64 %128, 1
  br i1 %125, label %357, label %127

127:                                              ; preds = %116, %122
  %128 = phi i64 [ 0, %116 ], [ %126, %122 ]
  %129 = phi i64 [ 2, %116 ], [ %274, %122 ]
  %130 = phi i64 [ 1, %116 ], [ %123, %122 ]
  %131 = phi i64 [ 3, %116 ], [ %124, %122 ]
  %132 = phi i32 [ 2, %116 ], [ %275, %122 ]
  %133 = add i64 %128, 1
  %134 = getelementptr double, double* %8, i64 %133
  %135 = getelementptr double, double* %9, i64 %133
  %136 = getelementptr double, double* %10, i64 %133
  %137 = getelementptr double, double* %11, i64 %133
  %138 = add i64 %128, 1
  %139 = add nsw i64 %129, -1
  %140 = getelementptr inbounds double, double* %8, i64 %139
  %141 = load double, double* %140, align 8, !tbaa !4
  %142 = getelementptr inbounds double, double* %9, i64 %139
  store double %141, double* %142, align 8, !tbaa !4
  %143 = trunc i64 %139 to i32
  %144 = shl i32 %143, 1
  %145 = or i32 %144, 1
  %146 = sitofp i32 %145 to double
  %147 = insertelement <2 x double> %118, double %141, i32 1
  %148 = insertelement <2 x double> undef, double %146, i32 0
  %149 = shufflevector <2 x double> %148, <2 x double> %49, <2 x i32> <i32 0, i32 2>
  %150 = fmul fast <2 x double> %147, %149
  %151 = insertelement <2 x double> undef, double %141, i32 0
  %152 = insertelement <2 x double> %151, double %146, i32 1
  %153 = fmul fast <2 x double> %150, %152
  %154 = bitcast double* %140 to <2 x double>*
  store <2 x double> %153, <2 x double>* %154, align 8, !tbaa !4
  %155 = getelementptr inbounds double, double* %10, i64 %139
  %156 = load double, double* %155, align 8, !tbaa !4
  %157 = getelementptr inbounds double, double* %11, i64 %139
  store double %156, double* %157, align 8, !tbaa !4
  %158 = fmul fast double %156, %45
  %159 = extractelement <2 x double> %150, i32 1
  %160 = fadd fast double %158, %159
  %161 = fmul fast double %160, %146
  store double %161, double* %155, align 8, !tbaa !4
  %162 = fmul fast double %141, %91
  %163 = fmul fast double %156, %46
  %164 = fsub fast double %162, %163
  %165 = fmul fast double %164, %146
  %166 = getelementptr inbounds double, double* %10, i64 %129
  store double %165, double* %166, align 8, !tbaa !4
  %167 = icmp ult i64 %138, 2
  br i1 %167, label %241, label %168

168:                                              ; preds = %127
  %169 = icmp ugt double* %135, %8
  %170 = icmp ugt double* %134, %9
  %171 = and i1 %169, %170
  %172 = icmp ugt double* %136, %8
  %173 = icmp ugt double* %134, %10
  %174 = and i1 %172, %173
  %175 = or i1 %171, %174
  %176 = icmp ugt double* %137, %8
  %177 = icmp ugt double* %134, %11
  %178 = and i1 %176, %177
  %179 = or i1 %175, %178
  %180 = icmp ugt double* %136, %9
  %181 = icmp ugt double* %135, %10
  %182 = and i1 %180, %181
  %183 = or i1 %179, %182
  %184 = icmp ugt double* %137, %9
  %185 = icmp ugt double* %135, %11
  %186 = and i1 %184, %185
  %187 = or i1 %183, %186
  %188 = icmp ugt double* %137, %10
  %189 = icmp ugt double* %136, %11
  %190 = and i1 %188, %189
  %191 = or i1 %187, %190
  br i1 %191, label %241, label %192

192:                                              ; preds = %168
  %193 = and i64 %138, -2
  %194 = shufflevector <2 x double> %150, <2 x double> undef, <2 x i32> zeroinitializer
  %195 = insertelement <2 x i64> undef, i64 %139, i32 0
  %196 = shufflevector <2 x i64> %195, <2 x i64> undef, <2 x i32> zeroinitializer
  %197 = insertelement <2 x i64> undef, i64 %129, i32 0
  %198 = shufflevector <2 x i64> %197, <2 x i64> undef, <2 x i32> zeroinitializer
  %199 = shufflevector <2 x double> %148, <2 x double> undef, <2 x i32> zeroinitializer
  br label %200

200:                                              ; preds = %200, %192
  %201 = phi i64 [ 0, %192 ], [ %236, %200 ]
  %202 = phi <2 x i64> [ <i64 0, i64 1>, %192 ], [ %237, %200 ]
  %203 = getelementptr inbounds double, double* %8, i64 %201
  %204 = bitcast double* %203 to <2 x double>*
  %205 = load <2 x double>, <2 x double>* %204, align 8, !tbaa !4, !alias.scope !18, !noalias !21
  %206 = fmul fast <2 x double> %194, %205
  %207 = add nuw nsw <2 x i64> %196, %202
  %208 = trunc <2 x i64> %207 to <2 x i32>
  %209 = sitofp <2 x i32> %208 to <2 x double>
  %210 = getelementptr inbounds double, double* %9, i64 %201
  %211 = bitcast double* %210 to <2 x double>*
  %212 = load <2 x double>, <2 x double>* %211, align 8, !tbaa !4, !alias.scope !25, !noalias !26
  %213 = fmul fast <2 x double> %212, %209
  %214 = fsub fast <2 x double> %206, %213
  %215 = sub nsw <2 x i64> %198, %202
  %216 = trunc <2 x i64> %215 to <2 x i32>
  %217 = sitofp <2 x i32> %216 to <2 x double>
  %218 = fdiv fast <2 x double> %214, %217
  %219 = bitcast double* %203 to <2 x double>*
  store <2 x double> %218, <2 x double>* %219, align 8, !tbaa !4, !alias.scope !18, !noalias !21
  %220 = bitcast double* %210 to <2 x double>*
  store <2 x double> %205, <2 x double>* %220, align 8, !tbaa !4, !alias.scope !25, !noalias !26
  %221 = getelementptr inbounds double, double* %10, i64 %201
  %222 = bitcast double* %221 to <2 x double>*
  %223 = load <2 x double>, <2 x double>* %222, align 8, !tbaa !4, !alias.scope !27, !noalias !28
  %224 = fmul fast <2 x double> %223, %120
  %225 = fmul fast <2 x double> %205, %121
  %226 = fsub fast <2 x double> %224, %225
  %227 = fmul fast <2 x double> %226, %199
  %228 = getelementptr inbounds double, double* %11, i64 %201
  %229 = bitcast double* %228 to <2 x double>*
  %230 = load <2 x double>, <2 x double>* %229, align 8, !tbaa !4, !alias.scope !28
  %231 = fmul fast <2 x double> %230, %209
  %232 = fsub fast <2 x double> %227, %231
  %233 = fdiv fast <2 x double> %232, %217
  %234 = bitcast double* %221 to <2 x double>*
  store <2 x double> %233, <2 x double>* %234, align 8, !tbaa !4, !alias.scope !27, !noalias !28
  %235 = bitcast double* %228 to <2 x double>*
  store <2 x double> %223, <2 x double>* %235, align 8, !tbaa !4, !alias.scope !28
  %236 = add i64 %201, 2
  %237 = add <2 x i64> %202, <i64 2, i64 2>
  %238 = icmp eq i64 %236, %193
  br i1 %238, label %239, label %200, !llvm.loop !29

239:                                              ; preds = %200
  %240 = icmp eq i64 %138, %193
  br i1 %240, label %273, label %241

241:                                              ; preds = %239, %168, %127
  %242 = phi i64 [ 0, %168 ], [ 0, %127 ], [ %193, %239 ]
  %243 = extractelement <2 x double> %150, i32 0
  br label %244

244:                                              ; preds = %241, %244
  %245 = phi i64 [ %271, %244 ], [ %242, %241 ]
  %246 = getelementptr inbounds double, double* %8, i64 %245
  %247 = load double, double* %246, align 8, !tbaa !4
  %248 = fmul fast double %243, %247
  %249 = add nuw nsw i64 %139, %245
  %250 = trunc i64 %249 to i32
  %251 = sitofp i32 %250 to double
  %252 = getelementptr inbounds double, double* %9, i64 %245
  %253 = load double, double* %252, align 8, !tbaa !4
  %254 = fmul fast double %253, %251
  %255 = fsub fast double %248, %254
  %256 = sub nsw i64 %129, %245
  %257 = trunc i64 %256 to i32
  %258 = sitofp i32 %257 to double
  %259 = fdiv fast double %255, %258
  store double %259, double* %246, align 8, !tbaa !4
  store double %247, double* %252, align 8, !tbaa !4
  %260 = getelementptr inbounds double, double* %10, i64 %245
  %261 = load double, double* %260, align 8, !tbaa !4
  %262 = fmul fast double %261, %45
  %263 = fmul fast double %247, %46
  %264 = fsub fast double %262, %263
  %265 = fmul fast double %264, %146
  %266 = getelementptr inbounds double, double* %11, i64 %245
  %267 = load double, double* %266, align 8, !tbaa !4
  %268 = fmul fast double %267, %251
  %269 = fsub fast double %265, %268
  %270 = fdiv fast double %269, %258
  store double %270, double* %260, align 8, !tbaa !4
  store double %261, double* %266, align 8, !tbaa !4
  %271 = add nuw nsw i64 %245, 1
  %272 = icmp eq i64 %271, %130
  br i1 %272, label %273, label %244, !llvm.loop !30

273:                                              ; preds = %244, %239
  %274 = add nuw nsw i64 %129, 1
  %275 = add nuw nsw i32 %132, 1
  %276 = trunc i64 %129 to i32
  %277 = mul nsw i32 %275, %276
  %278 = lshr i32 %277, 1
  %279 = shl i32 %276, 1
  %280 = or i32 %279, 1
  %281 = sitofp i32 %280 to double
  br label %282

282:                                              ; preds = %282, %273
  %283 = phi i64 [ %354, %282 ], [ 0, %273 ]
  %284 = phi i32 [ %355, %282 ], [ 0, %273 ]
  %285 = add nuw nsw i32 %284, %278
  %286 = mul nsw i32 %285, %16
  %287 = add nsw i32 %286, %117
  %288 = sub nsw i64 %129, %283
  %289 = getelementptr inbounds double, double* %12, i64 %288
  %290 = load double, double* %289, align 8, !tbaa !4
  %291 = fmul fast double %290, %281
  %292 = add nuw nsw i64 %283, %129
  %293 = getelementptr inbounds double, double* %12, i64 %292
  %294 = load double, double* %293, align 8, !tbaa !4
  %295 = fmul fast double %294, %20
  %296 = fdiv fast double %291, %295
  %297 = tail call fast double @llvm.sqrt.f64(double %296)
  %298 = getelementptr inbounds double, double* %13, i64 %283
  store double %297, double* %298, align 8, !tbaa !4
  %299 = trunc i64 %283 to i32
  %300 = sitofp i32 %299 to double
  %301 = load double, double* %73, align 8, !tbaa !4
  %302 = fmul fast double %301, %300
  %303 = tail call fast double @llvm.cos.f64(double %302)
  %304 = getelementptr inbounds double, double* %8, i64 %283
  %305 = load double, double* %304, align 8, !tbaa !4
  %306 = fmul fast double %297, %305
  %307 = fmul fast double %306, %303
  %308 = sext i32 %287 to i64
  %309 = getelementptr inbounds double, double* %0, i64 %308
  store double %307, double* %309, align 8, !tbaa !4
  %310 = load double, double* %298, align 8, !tbaa !4
  %311 = load double, double* %73, align 8, !tbaa !4
  %312 = fmul fast double %311, %300
  %313 = tail call fast double @llvm.sin.f64(double %312)
  %314 = load double, double* %304, align 8, !tbaa !4
  %315 = fmul fast double %314, %310
  %316 = fmul fast double %315, %313
  %317 = getelementptr inbounds double, double* %1, i64 %308
  store double %316, double* %317, align 8, !tbaa !4
  %318 = load double, double* %298, align 8, !tbaa !4
  %319 = load double, double* %73, align 8, !tbaa !4
  %320 = fmul fast double %319, %300
  %321 = tail call fast double @llvm.cos.f64(double %320)
  %322 = getelementptr inbounds double, double* %10, i64 %283
  %323 = load double, double* %322, align 8, !tbaa !4
  %324 = fmul fast double %323, %318
  %325 = fmul fast double %324, %321
  %326 = getelementptr inbounds double, double* %2, i64 %308
  store double %325, double* %326, align 8, !tbaa !4
  %327 = load double, double* %298, align 8, !tbaa !4
  %328 = load double, double* %73, align 8, !tbaa !4
  %329 = fmul fast double %328, %300
  %330 = tail call fast double @llvm.sin.f64(double %329)
  %331 = load double, double* %322, align 8, !tbaa !4
  %332 = fmul fast double %331, %327
  %333 = fmul fast double %332, %330
  %334 = getelementptr inbounds double, double* %3, i64 %308
  store double %333, double* %334, align 8, !tbaa !4
  %335 = load double, double* %298, align 8, !tbaa !4
  %336 = load double, double* %73, align 8, !tbaa !4
  %337 = fmul fast double %336, %300
  %338 = tail call fast double @llvm.sin.f64(double %337)
  %339 = load double, double* %304, align 8, !tbaa !4
  %340 = fsub fast double -0.000000e+00, %300
  %341 = fmul fast double %335, %340
  %342 = fmul fast double %341, %339
  %343 = fmul fast double %342, %338
  %344 = getelementptr inbounds double, double* %4, i64 %308
  store double %343, double* %344, align 8, !tbaa !4
  %345 = load double, double* %298, align 8, !tbaa !4
  %346 = load double, double* %73, align 8, !tbaa !4
  %347 = fmul fast double %346, %300
  %348 = tail call fast double @llvm.cos.f64(double %347)
  %349 = load double, double* %304, align 8, !tbaa !4
  %350 = fmul fast double %345, %300
  %351 = fmul fast double %350, %349
  %352 = fmul fast double %351, %348
  %353 = getelementptr inbounds double, double* %5, i64 %308
  store double %352, double* %353, align 8, !tbaa !4
  %354 = add nuw nsw i64 %283, 1
  %355 = add nuw nsw i32 %284, 1
  %356 = icmp eq i64 %354, %131
  br i1 %356, label %122, label %282

357:                                              ; preds = %122, %35
  %358 = add nuw nsw i64 %36, 1
  %359 = icmp eq i64 %358, %31
  br i1 %359, label %34, label %35
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z26cpuSphericalHarmonicsDerivIfEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S0_ii(float* %0, float* %1, float* %2, float* %3, float* %4, float* %5, float* %6, float* %7, float* %8, float* %9, float* %10, float* %11, float* %12, float* %13, float %14, i32 %15, i32 %16) local_unnamed_addr #0 {
  %18 = icmp sgt i32 %16, 0
  br i1 %18, label %19, label %34

19:                                               ; preds = %17
  %20 = fmul fast float %14, 4.000000e+00
  %21 = fdiv fast float 2.500000e-01, %14
  %22 = tail call fast float @llvm.sqrt.f32(float %21) #8
  %23 = getelementptr inbounds float, float* %8, i64 1
  %24 = getelementptr inbounds float, float* %10, i64 1
  %25 = getelementptr inbounds float, float* %12, i64 1
  %26 = getelementptr inbounds float, float* %12, i64 2
  %27 = getelementptr inbounds float, float* %13, i64 1
  %28 = shl nuw i32 %16, 1
  %29 = icmp slt i32 %15, 2
  %30 = add i32 %15, 1
  %31 = sext i32 %28 to i64
  %32 = zext i32 %16 to i64
  %33 = zext i32 %30 to i64
  br label %35

34:                                               ; preds = %351, %17
  ret void

35:                                               ; preds = %351, %19
  %36 = phi i64 [ 0, %19 ], [ %352, %351 ]
  %37 = getelementptr inbounds float, float* %0, i64 %36
  store float %22, float* %37, align 4, !tbaa !8
  %38 = getelementptr inbounds float, float* %1, i64 %36
  store float 0.000000e+00, float* %38, align 4, !tbaa !8
  %39 = getelementptr inbounds float, float* %2, i64 %36
  store float 0.000000e+00, float* %39, align 4, !tbaa !8
  %40 = getelementptr inbounds float, float* %4, i64 %36
  store float 0.000000e+00, float* %40, align 4, !tbaa !8
  %41 = getelementptr inbounds float, float* %3, i64 %36
  store float 0.000000e+00, float* %41, align 4, !tbaa !8
  %42 = getelementptr inbounds float, float* %5, i64 %36
  store float 0.000000e+00, float* %42, align 4, !tbaa !8
  %43 = getelementptr inbounds float, float* %6, i64 %36
  %44 = load float, float* %43, align 4, !tbaa !8
  %45 = tail call fast float @llvm.cos.f32(float %44) #8
  %46 = tail call fast float @llvm.sin.f32(float %44) #8
  %47 = fneg fast float %46
  %48 = fneg fast float %45
  store float %45, float* %8, align 4, !tbaa !8
  store float %47, float* %23, align 4, !tbaa !8
  store float %47, float* %10, align 4, !tbaa !8
  store float %48, float* %24, align 4, !tbaa !8
  %49 = load float, float* %25, align 4, !tbaa !8
  %50 = fmul fast float %49, 3.000000e+00
  %51 = fmul fast float %49, %20
  %52 = fdiv fast float %50, %51
  %53 = tail call fast float @llvm.sqrt.f32(float %52) #8
  store float %53, float* %13, align 4, !tbaa !8
  %54 = load float, float* %8, align 4, !tbaa !8
  %55 = fmul fast float %53, %54
  %56 = add nuw nsw i64 %36, %32
  %57 = getelementptr inbounds float, float* %0, i64 %56
  store float %55, float* %57, align 4, !tbaa !8
  %58 = getelementptr inbounds float, float* %1, i64 %56
  store float 0.000000e+00, float* %58, align 4, !tbaa !8
  %59 = load float, float* %13, align 4, !tbaa !8
  %60 = fmul fast float %59, %47
  %61 = getelementptr inbounds float, float* %2, i64 %56
  store float %60, float* %61, align 4, !tbaa !8
  %62 = getelementptr inbounds float, float* %3, i64 %56
  store float 0.000000e+00, float* %62, align 4, !tbaa !8
  %63 = getelementptr inbounds float, float* %4, i64 %56
  store float 0.000000e+00, float* %63, align 4, !tbaa !8
  %64 = getelementptr inbounds float, float* %5, i64 %56
  store float 0.000000e+00, float* %64, align 4, !tbaa !8
  %65 = load float, float* %12, align 4, !tbaa !8
  %66 = fmul fast float %65, 3.000000e+00
  %67 = load float, float* %26, align 4, !tbaa !8
  %68 = fmul fast float %67, %20
  %69 = fdiv fast float %66, %68
  %70 = tail call fast float @llvm.sqrt.f32(float %69) #8
  store float %70, float* %27, align 4, !tbaa !8
  %71 = getelementptr inbounds float, float* %7, i64 %36
  %72 = load float, float* %71, align 4, !tbaa !8
  %73 = tail call fast float @llvm.cos.f32(float %72) #8
  %74 = fmul fast float %70, %73
  %75 = load float, float* %23, align 4, !tbaa !8
  %76 = fmul fast float %74, %75
  %77 = add nsw i64 %36, %31
  %78 = getelementptr inbounds float, float* %0, i64 %77
  store float %76, float* %78, align 4, !tbaa !8
  %79 = load float, float* %27, align 4, !tbaa !8
  %80 = load float, float* %71, align 4, !tbaa !8
  %81 = tail call fast float @llvm.sin.f32(float %80) #8
  %82 = fmul fast float %81, %79
  %83 = load float, float* %23, align 4, !tbaa !8
  %84 = fmul fast float %82, %83
  %85 = getelementptr inbounds float, float* %1, i64 %77
  store float %84, float* %85, align 4, !tbaa !8
  %86 = load float, float* %27, align 4, !tbaa !8
  %87 = load float, float* %71, align 4, !tbaa !8
  %88 = tail call fast float @llvm.cos.f32(float %87) #8
  %89 = fmul fast float %86, %48
  %90 = fmul fast float %89, %88
  %91 = getelementptr inbounds float, float* %2, i64 %77
  store float %90, float* %91, align 4, !tbaa !8
  %92 = load float, float* %27, align 4, !tbaa !8
  %93 = load float, float* %71, align 4, !tbaa !8
  %94 = tail call fast float @llvm.sin.f32(float %93) #8
  %95 = fmul fast float %92, %48
  %96 = fmul fast float %95, %94
  %97 = getelementptr inbounds float, float* %3, i64 %77
  store float %96, float* %97, align 4, !tbaa !8
  %98 = load float, float* %27, align 4, !tbaa !8
  %99 = load float, float* %71, align 4, !tbaa !8
  %100 = tail call fast float @llvm.sin.f32(float %99) #8
  %101 = load float, float* %23, align 4, !tbaa !8
  %102 = fsub fast float -0.000000e+00, %98
  %103 = fmul fast float %101, %102
  %104 = fmul fast float %103, %100
  %105 = getelementptr inbounds float, float* %4, i64 %77
  store float %104, float* %105, align 4, !tbaa !8
  %106 = load float, float* %27, align 4, !tbaa !8
  %107 = load float, float* %71, align 4, !tbaa !8
  %108 = tail call fast float @llvm.cos.f32(float %107) #8
  %109 = load float, float* %23, align 4, !tbaa !8
  %110 = fmul fast float %109, %106
  %111 = fmul fast float %110, %108
  %112 = getelementptr inbounds float, float* %5, i64 %77
  store float %111, float* %112, align 4, !tbaa !8
  br i1 %29, label %351, label %113

113:                                              ; preds = %35
  %114 = trunc i64 %36 to i32
  %115 = insertelement <4 x float> undef, float %45, i32 0
  %116 = shufflevector <4 x float> %115, <4 x float> undef, <4 x i32> zeroinitializer
  %117 = insertelement <4 x float> undef, float %46, i32 0
  %118 = shufflevector <4 x float> %117, <4 x float> undef, <4 x i32> zeroinitializer
  br label %124

119:                                              ; preds = %276
  %120 = add nuw nsw i64 %127, 1
  %121 = add nuw nsw i64 %128, 1
  %122 = icmp eq i64 %268, %33
  %123 = add i64 %125, 1
  br i1 %122, label %351, label %124

124:                                              ; preds = %113, %119
  %125 = phi i64 [ 0, %113 ], [ %123, %119 ]
  %126 = phi i64 [ 2, %113 ], [ %268, %119 ]
  %127 = phi i64 [ 1, %113 ], [ %120, %119 ]
  %128 = phi i64 [ 3, %113 ], [ %121, %119 ]
  %129 = phi i32 [ 2, %113 ], [ %269, %119 ]
  %130 = add i64 %125, 1
  %131 = getelementptr float, float* %8, i64 %130
  %132 = getelementptr float, float* %9, i64 %130
  %133 = getelementptr float, float* %10, i64 %130
  %134 = getelementptr float, float* %11, i64 %130
  %135 = add i64 %125, 1
  %136 = add nsw i64 %126, -1
  %137 = getelementptr inbounds float, float* %8, i64 %136
  %138 = load float, float* %137, align 4, !tbaa !8
  %139 = getelementptr inbounds float, float* %9, i64 %136
  store float %138, float* %139, align 4, !tbaa !8
  %140 = trunc i64 %136 to i32
  %141 = shl i32 %140, 1
  %142 = or i32 %141, 1
  %143 = sitofp i32 %142 to float
  %144 = fmul fast float %45, %143
  %145 = fmul fast float %144, %138
  store float %145, float* %137, align 4, !tbaa !8
  %146 = fmul fast float %138, %47
  %147 = fmul fast float %146, %143
  %148 = getelementptr inbounds float, float* %8, i64 %126
  store float %147, float* %148, align 4, !tbaa !8
  %149 = getelementptr inbounds float, float* %10, i64 %136
  %150 = load float, float* %149, align 4, !tbaa !8
  %151 = getelementptr inbounds float, float* %11, i64 %136
  store float %150, float* %151, align 4, !tbaa !8
  %152 = fmul fast float %150, %45
  %153 = fadd fast float %152, %146
  %154 = fmul fast float %153, %143
  store float %154, float* %149, align 4, !tbaa !8
  %155 = fmul fast float %138, %48
  %156 = fmul fast float %150, %46
  %157 = fsub fast float %155, %156
  %158 = fmul fast float %157, %143
  %159 = getelementptr inbounds float, float* %10, i64 %126
  store float %158, float* %159, align 4, !tbaa !8
  %160 = icmp ult i64 %135, 4
  br i1 %160, label %236, label %161

161:                                              ; preds = %124
  %162 = icmp ugt float* %132, %8
  %163 = icmp ugt float* %131, %9
  %164 = and i1 %162, %163
  %165 = icmp ugt float* %133, %8
  %166 = icmp ugt float* %131, %10
  %167 = and i1 %165, %166
  %168 = or i1 %164, %167
  %169 = icmp ugt float* %134, %8
  %170 = icmp ugt float* %131, %11
  %171 = and i1 %169, %170
  %172 = or i1 %168, %171
  %173 = icmp ugt float* %133, %9
  %174 = icmp ugt float* %132, %10
  %175 = and i1 %173, %174
  %176 = or i1 %172, %175
  %177 = icmp ugt float* %134, %9
  %178 = icmp ugt float* %132, %11
  %179 = and i1 %177, %178
  %180 = or i1 %176, %179
  %181 = icmp ugt float* %134, %10
  %182 = icmp ugt float* %133, %11
  %183 = and i1 %181, %182
  %184 = or i1 %180, %183
  br i1 %184, label %236, label %185

185:                                              ; preds = %161
  %186 = and i64 %135, -4
  %187 = insertelement <4 x float> undef, float %144, i32 0
  %188 = shufflevector <4 x float> %187, <4 x float> undef, <4 x i32> zeroinitializer
  %189 = insertelement <4 x i64> undef, i64 %136, i32 0
  %190 = shufflevector <4 x i64> %189, <4 x i64> undef, <4 x i32> zeroinitializer
  %191 = insertelement <4 x i64> undef, i64 %126, i32 0
  %192 = shufflevector <4 x i64> %191, <4 x i64> undef, <4 x i32> zeroinitializer
  %193 = insertelement <4 x float> undef, float %143, i32 0
  %194 = shufflevector <4 x float> %193, <4 x float> undef, <4 x i32> zeroinitializer
  br label %195

195:                                              ; preds = %195, %185
  %196 = phi i64 [ 0, %185 ], [ %231, %195 ]
  %197 = phi <4 x i64> [ <i64 0, i64 1, i64 2, i64 3>, %185 ], [ %232, %195 ]
  %198 = getelementptr inbounds float, float* %8, i64 %196
  %199 = bitcast float* %198 to <4 x float>*
  %200 = load <4 x float>, <4 x float>* %199, align 4, !tbaa !8, !alias.scope !31, !noalias !34
  %201 = fmul fast <4 x float> %188, %200
  %202 = add nuw nsw <4 x i64> %190, %197
  %203 = trunc <4 x i64> %202 to <4 x i32>
  %204 = sitofp <4 x i32> %203 to <4 x float>
  %205 = getelementptr inbounds float, float* %9, i64 %196
  %206 = bitcast float* %205 to <4 x float>*
  %207 = load <4 x float>, <4 x float>* %206, align 4, !tbaa !8, !alias.scope !38, !noalias !39
  %208 = fmul fast <4 x float> %207, %204
  %209 = fsub fast <4 x float> %201, %208
  %210 = sub nsw <4 x i64> %192, %197
  %211 = trunc <4 x i64> %210 to <4 x i32>
  %212 = sitofp <4 x i32> %211 to <4 x float>
  %213 = fdiv fast <4 x float> %209, %212
  %214 = bitcast float* %198 to <4 x float>*
  store <4 x float> %213, <4 x float>* %214, align 4, !tbaa !8, !alias.scope !31, !noalias !34
  %215 = bitcast float* %205 to <4 x float>*
  store <4 x float> %200, <4 x float>* %215, align 4, !tbaa !8, !alias.scope !38, !noalias !39
  %216 = getelementptr inbounds float, float* %10, i64 %196
  %217 = bitcast float* %216 to <4 x float>*
  %218 = load <4 x float>, <4 x float>* %217, align 4, !tbaa !8, !alias.scope !40, !noalias !41
  %219 = fmul fast <4 x float> %218, %116
  %220 = fmul fast <4 x float> %200, %118
  %221 = fsub fast <4 x float> %219, %220
  %222 = fmul fast <4 x float> %221, %194
  %223 = getelementptr inbounds float, float* %11, i64 %196
  %224 = bitcast float* %223 to <4 x float>*
  %225 = load <4 x float>, <4 x float>* %224, align 4, !tbaa !8, !alias.scope !41
  %226 = fmul fast <4 x float> %225, %204
  %227 = fsub fast <4 x float> %222, %226
  %228 = fdiv fast <4 x float> %227, %212
  %229 = bitcast float* %216 to <4 x float>*
  store <4 x float> %228, <4 x float>* %229, align 4, !tbaa !8, !alias.scope !40, !noalias !41
  %230 = bitcast float* %223 to <4 x float>*
  store <4 x float> %218, <4 x float>* %230, align 4, !tbaa !8, !alias.scope !41
  %231 = add i64 %196, 4
  %232 = add <4 x i64> %197, <i64 4, i64 4, i64 4, i64 4>
  %233 = icmp eq i64 %231, %186
  br i1 %233, label %234, label %195, !llvm.loop !42

234:                                              ; preds = %195
  %235 = icmp eq i64 %135, %186
  br i1 %235, label %267, label %236

236:                                              ; preds = %234, %161, %124
  %237 = phi i64 [ 0, %161 ], [ 0, %124 ], [ %186, %234 ]
  br label %238

238:                                              ; preds = %236, %238
  %239 = phi i64 [ %265, %238 ], [ %237, %236 ]
  %240 = getelementptr inbounds float, float* %8, i64 %239
  %241 = load float, float* %240, align 4, !tbaa !8
  %242 = fmul fast float %144, %241
  %243 = add nuw nsw i64 %136, %239
  %244 = trunc i64 %243 to i32
  %245 = sitofp i32 %244 to float
  %246 = getelementptr inbounds float, float* %9, i64 %239
  %247 = load float, float* %246, align 4, !tbaa !8
  %248 = fmul fast float %247, %245
  %249 = fsub fast float %242, %248
  %250 = sub nsw i64 %126, %239
  %251 = trunc i64 %250 to i32
  %252 = sitofp i32 %251 to float
  %253 = fdiv fast float %249, %252
  store float %253, float* %240, align 4, !tbaa !8
  store float %241, float* %246, align 4, !tbaa !8
  %254 = getelementptr inbounds float, float* %10, i64 %239
  %255 = load float, float* %254, align 4, !tbaa !8
  %256 = fmul fast float %255, %45
  %257 = fmul fast float %241, %46
  %258 = fsub fast float %256, %257
  %259 = fmul fast float %258, %143
  %260 = getelementptr inbounds float, float* %11, i64 %239
  %261 = load float, float* %260, align 4, !tbaa !8
  %262 = fmul fast float %261, %245
  %263 = fsub fast float %259, %262
  %264 = fdiv fast float %263, %252
  store float %264, float* %254, align 4, !tbaa !8
  store float %255, float* %260, align 4, !tbaa !8
  %265 = add nuw nsw i64 %239, 1
  %266 = icmp eq i64 %265, %127
  br i1 %266, label %267, label %238, !llvm.loop !43

267:                                              ; preds = %238, %234
  %268 = add nuw nsw i64 %126, 1
  %269 = add nuw nsw i32 %129, 1
  %270 = trunc i64 %126 to i32
  %271 = mul nsw i32 %269, %270
  %272 = lshr i32 %271, 1
  %273 = shl i32 %270, 1
  %274 = or i32 %273, 1
  %275 = sitofp i32 %274 to float
  br label %276

276:                                              ; preds = %276, %267
  %277 = phi i64 [ %348, %276 ], [ 0, %267 ]
  %278 = phi i32 [ %349, %276 ], [ 0, %267 ]
  %279 = add nuw nsw i32 %278, %272
  %280 = mul nsw i32 %279, %16
  %281 = add nsw i32 %280, %114
  %282 = sub nsw i64 %126, %277
  %283 = getelementptr inbounds float, float* %12, i64 %282
  %284 = load float, float* %283, align 4, !tbaa !8
  %285 = fmul fast float %284, %275
  %286 = add nuw nsw i64 %277, %126
  %287 = getelementptr inbounds float, float* %12, i64 %286
  %288 = load float, float* %287, align 4, !tbaa !8
  %289 = fmul fast float %288, %20
  %290 = fdiv fast float %285, %289
  %291 = tail call fast float @llvm.sqrt.f32(float %290) #8
  %292 = getelementptr inbounds float, float* %13, i64 %277
  store float %291, float* %292, align 4, !tbaa !8
  %293 = trunc i64 %277 to i32
  %294 = sitofp i32 %293 to float
  %295 = load float, float* %71, align 4, !tbaa !8
  %296 = fmul fast float %295, %294
  %297 = tail call fast float @llvm.cos.f32(float %296) #8
  %298 = getelementptr inbounds float, float* %8, i64 %277
  %299 = load float, float* %298, align 4, !tbaa !8
  %300 = fmul fast float %291, %299
  %301 = fmul fast float %300, %297
  %302 = sext i32 %281 to i64
  %303 = getelementptr inbounds float, float* %0, i64 %302
  store float %301, float* %303, align 4, !tbaa !8
  %304 = load float, float* %292, align 4, !tbaa !8
  %305 = load float, float* %71, align 4, !tbaa !8
  %306 = fmul fast float %305, %294
  %307 = tail call fast float @llvm.sin.f32(float %306) #8
  %308 = load float, float* %298, align 4, !tbaa !8
  %309 = fmul fast float %308, %304
  %310 = fmul fast float %309, %307
  %311 = getelementptr inbounds float, float* %1, i64 %302
  store float %310, float* %311, align 4, !tbaa !8
  %312 = load float, float* %292, align 4, !tbaa !8
  %313 = load float, float* %71, align 4, !tbaa !8
  %314 = fmul fast float %313, %294
  %315 = tail call fast float @llvm.cos.f32(float %314) #8
  %316 = getelementptr inbounds float, float* %10, i64 %277
  %317 = load float, float* %316, align 4, !tbaa !8
  %318 = fmul fast float %317, %312
  %319 = fmul fast float %318, %315
  %320 = getelementptr inbounds float, float* %2, i64 %302
  store float %319, float* %320, align 4, !tbaa !8
  %321 = load float, float* %292, align 4, !tbaa !8
  %322 = load float, float* %71, align 4, !tbaa !8
  %323 = fmul fast float %322, %294
  %324 = tail call fast float @llvm.sin.f32(float %323) #8
  %325 = load float, float* %316, align 4, !tbaa !8
  %326 = fmul fast float %325, %321
  %327 = fmul fast float %326, %324
  %328 = getelementptr inbounds float, float* %3, i64 %302
  store float %327, float* %328, align 4, !tbaa !8
  %329 = load float, float* %292, align 4, !tbaa !8
  %330 = load float, float* %71, align 4, !tbaa !8
  %331 = fmul fast float %330, %294
  %332 = tail call fast float @llvm.sin.f32(float %331) #8
  %333 = load float, float* %298, align 4, !tbaa !8
  %334 = fsub fast float -0.000000e+00, %294
  %335 = fmul fast float %329, %334
  %336 = fmul fast float %335, %333
  %337 = fmul fast float %336, %332
  %338 = getelementptr inbounds float, float* %4, i64 %302
  store float %337, float* %338, align 4, !tbaa !8
  %339 = load float, float* %292, align 4, !tbaa !8
  %340 = load float, float* %71, align 4, !tbaa !8
  %341 = fmul fast float %340, %294
  %342 = tail call fast float @llvm.cos.f32(float %341) #8
  %343 = load float, float* %298, align 4, !tbaa !8
  %344 = fmul fast float %339, %294
  %345 = fmul fast float %344, %343
  %346 = fmul fast float %345, %342
  %347 = getelementptr inbounds float, float* %5, i64 %302
  store float %346, float* %347, align 4, !tbaa !8
  %348 = add nuw nsw i64 %277, 1
  %349 = add nuw nsw i32 %278, 1
  %350 = icmp eq i64 %348, %128
  br i1 %350, label %119, label %276

351:                                              ; preds = %119, %35
  %352 = add nuw nsw i64 %36, 1
  %353 = icmp eq i64 %352, %32
  br i1 %353, label %34, label %35
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z18cpuSphericalBesselIdEvPT_S1_S1_S1_iii(double* %0, double* %1, double* %2, double* %3, i32 %4, i32 %5, i32 %6) local_unnamed_addr #0 {
  %8 = icmp sgt i32 %6, 0
  br i1 %8, label %9, label %25

9:                                                ; preds = %7
  %10 = icmp sgt i32 %5, 0
  %11 = mul i32 %6, %5
  %12 = icmp slt i32 %4, 2
  %13 = getelementptr inbounds double, double* %3, i64 1
  %14 = sext i32 %5 to i64
  %15 = add i32 %4, 1
  %16 = zext i32 %6 to i64
  %17 = zext i32 %5 to i64
  %18 = zext i32 %15 to i64
  %19 = getelementptr inbounds double, double* %3, i64 2
  br label %20

20:                                               ; preds = %108, %9
  %21 = phi i64 [ 0, %9 ], [ %109, %108 ]
  br i1 %10, label %22, label %108

22:                                               ; preds = %20
  %23 = getelementptr inbounds double, double* %1, i64 %21
  %24 = trunc i64 %21 to i32
  br label %26

25:                                               ; preds = %108, %7
  ret void

26:                                               ; preds = %105, %22
  %27 = phi i64 [ 0, %22 ], [ %106, %105 ]
  %28 = getelementptr inbounds double, double* %2, i64 %27
  %29 = load double, double* %28, align 8, !tbaa !4
  %30 = load double, double* %23, align 8, !tbaa !4
  %31 = fmul fast double %30, %29
  %32 = tail call fast double @llvm.cos.f64(double %31)
  %33 = fdiv fast double %32, %31
  %34 = mul nsw i64 %27, %16
  %35 = add nuw nsw i64 %34, %21
  %36 = getelementptr inbounds double, double* %0, i64 %35
  store double %33, double* %36, align 8, !tbaa !4
  %37 = add nsw i64 %27, %14
  %38 = getelementptr inbounds double, double* %2, i64 %37
  %39 = load double, double* %38, align 8, !tbaa !4
  %40 = load double, double* %23, align 8, !tbaa !4
  %41 = fmul fast double %40, %39
  %42 = tail call fast double @llvm.cos.f64(double %41)
  %43 = fmul fast double %41, %41
  %44 = fdiv fast double %42, %43
  %45 = tail call fast double @llvm.sin.f64(double %41)
  %46 = fdiv fast double %45, %41
  %47 = fadd fast double %44, %46
  %48 = trunc i64 %34 to i32
  %49 = add i32 %48, %24
  %50 = add i32 %49, %11
  %51 = sext i32 %50 to i64
  %52 = getelementptr inbounds double, double* %0, i64 %51
  store double %47, double* %52, align 8, !tbaa !4
  br i1 %12, label %105, label %53

53:                                               ; preds = %26, %76
  %54 = phi i64 [ %77, %76 ], [ 2, %26 ]
  %55 = phi i64 [ %86, %76 ], [ 3, %26 ]
  %56 = mul nsw i64 %54, %14
  %57 = add nsw i64 %56, %27
  %58 = getelementptr inbounds double, double* %2, i64 %57
  %59 = load double, double* %58, align 8, !tbaa !4
  %60 = load double, double* %23, align 8, !tbaa !4
  %61 = fmul fast double %60, %59
  %62 = tail call fast double @llvm.cos.f64(double %61)
  %63 = fneg fast double %62
  %64 = fdiv fast double %63, %61
  store double %64, double* %3, align 8, !tbaa !4
  %65 = fmul fast double %61, %61
  %66 = fdiv fast double %63, %65
  %67 = tail call fast double @llvm.sin.f64(double %61)
  %68 = fdiv fast double %67, %61
  %69 = fsub fast double %66, %68
  store double %69, double* %13, align 8, !tbaa !4
  %70 = fmul fast double %69, 3.000000e+00
  %71 = fdiv fast double %70, %61
  %72 = fsub fast double %71, %64
  store double %72, double* %19, align 8, !tbaa !4
  %73 = icmp eq i64 %55, 3
  br i1 %73, label %76, label %74

74:                                               ; preds = %53
  %75 = fdiv fast double 1.000000e+00, %61
  br label %88

76:                                               ; preds = %88, %53
  %77 = add nuw nsw i64 %54, 1
  %78 = getelementptr inbounds double, double* %3, i64 %54
  %79 = load double, double* %78, align 8, !tbaa !4
  %80 = fneg fast double %79
  %81 = trunc i64 %54 to i32
  %82 = mul i32 %11, %81
  %83 = add i32 %49, %82
  %84 = sext i32 %83 to i64
  %85 = getelementptr inbounds double, double* %0, i64 %84
  store double %80, double* %85, align 8, !tbaa !4
  %86 = add nuw nsw i64 %55, 1
  %87 = icmp eq i64 %77, %18
  br i1 %87, label %105, label %53

88:                                               ; preds = %74, %88
  %89 = phi double [ %102, %88 ], [ %72, %74 ]
  %90 = phi i64 [ %95, %88 ], [ 3, %74 ]
  %91 = phi i64 [ %90, %88 ], [ 2, %74 ]
  %92 = add nsw i64 %91, -1
  %93 = getelementptr inbounds double, double* %3, i64 %92
  %94 = load double, double* %93, align 8, !tbaa !4
  %95 = add nuw nsw i64 %90, 1
  %96 = trunc i64 %95 to i32
  %97 = shl nuw nsw i32 %96, 1
  %98 = add nsw i32 %97, -3
  %99 = sitofp i32 %98 to double
  %100 = fmul fast double %89, %99
  %101 = fmul fast double %100, %75
  %102 = fsub fast double %101, %94
  %103 = getelementptr inbounds double, double* %3, i64 %90
  store double %102, double* %103, align 8, !tbaa !4
  %104 = icmp eq i64 %95, %55
  br i1 %104, label %76, label %88

105:                                              ; preds = %76, %26
  %106 = add nuw nsw i64 %27, 1
  %107 = icmp eq i64 %106, %17
  br i1 %107, label %108, label %26

108:                                              ; preds = %105, %20
  %109 = add nuw nsw i64 %21, 1
  %110 = icmp eq i64 %109, %16
  br i1 %110, label %25, label %20
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z18cpuSphericalBesselIfEvPT_S1_S1_S1_iii(float* %0, float* %1, float* %2, float* %3, i32 %4, i32 %5, i32 %6) local_unnamed_addr #0 {
  %8 = icmp sgt i32 %6, 0
  br i1 %8, label %9, label %25

9:                                                ; preds = %7
  %10 = icmp sgt i32 %5, 0
  %11 = mul i32 %6, %5
  %12 = icmp slt i32 %4, 2
  %13 = getelementptr inbounds float, float* %3, i64 1
  %14 = sext i32 %5 to i64
  %15 = add i32 %4, 1
  %16 = zext i32 %6 to i64
  %17 = zext i32 %5 to i64
  %18 = zext i32 %15 to i64
  %19 = getelementptr inbounds float, float* %3, i64 2
  br label %20

20:                                               ; preds = %108, %9
  %21 = phi i64 [ 0, %9 ], [ %109, %108 ]
  br i1 %10, label %22, label %108

22:                                               ; preds = %20
  %23 = getelementptr inbounds float, float* %1, i64 %21
  %24 = trunc i64 %21 to i32
  br label %26

25:                                               ; preds = %108, %7
  ret void

26:                                               ; preds = %105, %22
  %27 = phi i64 [ 0, %22 ], [ %106, %105 ]
  %28 = getelementptr inbounds float, float* %2, i64 %27
  %29 = load float, float* %28, align 4, !tbaa !8
  %30 = load float, float* %23, align 4, !tbaa !8
  %31 = fmul fast float %30, %29
  %32 = tail call fast float @llvm.cos.f32(float %31) #8
  %33 = fdiv fast float %32, %31
  %34 = mul nsw i64 %27, %16
  %35 = add nuw nsw i64 %34, %21
  %36 = getelementptr inbounds float, float* %0, i64 %35
  store float %33, float* %36, align 4, !tbaa !8
  %37 = add nsw i64 %27, %14
  %38 = getelementptr inbounds float, float* %2, i64 %37
  %39 = load float, float* %38, align 4, !tbaa !8
  %40 = load float, float* %23, align 4, !tbaa !8
  %41 = fmul fast float %40, %39
  %42 = tail call fast float @llvm.cos.f32(float %41) #8
  %43 = fmul fast float %41, %41
  %44 = fdiv fast float %42, %43
  %45 = tail call fast float @llvm.sin.f32(float %41) #8
  %46 = fdiv fast float %45, %41
  %47 = fadd fast float %44, %46
  %48 = trunc i64 %34 to i32
  %49 = add i32 %48, %24
  %50 = add i32 %49, %11
  %51 = sext i32 %50 to i64
  %52 = getelementptr inbounds float, float* %0, i64 %51
  store float %47, float* %52, align 4, !tbaa !8
  br i1 %12, label %105, label %53

53:                                               ; preds = %26, %76
  %54 = phi i64 [ %77, %76 ], [ 2, %26 ]
  %55 = phi i64 [ %86, %76 ], [ 3, %26 ]
  %56 = mul nsw i64 %54, %14
  %57 = add nsw i64 %56, %27
  %58 = getelementptr inbounds float, float* %2, i64 %57
  %59 = load float, float* %58, align 4, !tbaa !8
  %60 = load float, float* %23, align 4, !tbaa !8
  %61 = fmul fast float %60, %59
  %62 = tail call fast float @llvm.cos.f32(float %61) #8
  %63 = fneg fast float %62
  %64 = fdiv fast float %63, %61
  store float %64, float* %3, align 4, !tbaa !8
  %65 = fmul fast float %61, %61
  %66 = fdiv fast float %63, %65
  %67 = tail call fast float @llvm.sin.f32(float %61) #8
  %68 = fdiv fast float %67, %61
  %69 = fsub fast float %66, %68
  store float %69, float* %13, align 4, !tbaa !8
  %70 = fmul fast float %69, 3.000000e+00
  %71 = fdiv fast float %70, %61
  %72 = fsub fast float %71, %64
  store float %72, float* %19, align 4, !tbaa !8
  %73 = icmp eq i64 %55, 3
  br i1 %73, label %76, label %74

74:                                               ; preds = %53
  %75 = fdiv fast float 1.000000e+00, %61
  br label %88

76:                                               ; preds = %88, %53
  %77 = add nuw nsw i64 %54, 1
  %78 = getelementptr inbounds float, float* %3, i64 %54
  %79 = load float, float* %78, align 4, !tbaa !8
  %80 = fneg fast float %79
  %81 = trunc i64 %54 to i32
  %82 = mul i32 %11, %81
  %83 = add i32 %49, %82
  %84 = sext i32 %83 to i64
  %85 = getelementptr inbounds float, float* %0, i64 %84
  store float %80, float* %85, align 4, !tbaa !8
  %86 = add nuw nsw i64 %55, 1
  %87 = icmp eq i64 %77, %18
  br i1 %87, label %105, label %53

88:                                               ; preds = %74, %88
  %89 = phi float [ %102, %88 ], [ %72, %74 ]
  %90 = phi i64 [ %95, %88 ], [ 3, %74 ]
  %91 = phi i64 [ %90, %88 ], [ 2, %74 ]
  %92 = add nsw i64 %91, -1
  %93 = getelementptr inbounds float, float* %3, i64 %92
  %94 = load float, float* %93, align 4, !tbaa !8
  %95 = add nuw nsw i64 %90, 1
  %96 = trunc i64 %95 to i32
  %97 = shl nuw nsw i32 %96, 1
  %98 = add nsw i32 %97, -3
  %99 = sitofp i32 %98 to float
  %100 = fmul fast float %89, %99
  %101 = fmul fast float %100, %75
  %102 = fsub fast float %101, %94
  %103 = getelementptr inbounds float, float* %3, i64 %90
  store float %102, float* %103, align 4, !tbaa !8
  %104 = icmp eq i64 %95, %55
  br i1 %104, label %76, label %88

105:                                              ; preds = %76, %26
  %106 = add nuw nsw i64 %27, 1
  %107 = icmp eq i64 %106, %17
  br i1 %107, label %108, label %26

108:                                              ; preds = %105, %20
  %109 = add nuw nsw i64 %21, 1
  %110 = icmp eq i64 %109, %16
  br i1 %110, label %25, label %20
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z23cpuSphericalBesselDerivIdEvPT_S1_S1_S1_S1_S1_iii(double* %0, double* %1, double* %2, double* %3, double* %4, double* %5, i32 %6, i32 %7, i32 %8) local_unnamed_addr #0 {
  %10 = icmp sgt i32 %8, 0
  br i1 %10, label %11, label %29

11:                                               ; preds = %9
  %12 = icmp sgt i32 %7, 0
  %13 = mul i32 %8, %7
  %14 = icmp slt i32 %6, 2
  %15 = getelementptr inbounds double, double* %4, i64 1
  %16 = getelementptr inbounds double, double* %5, i64 1
  %17 = sext i32 %7 to i64
  %18 = add i32 %6, 1
  %19 = zext i32 %8 to i64
  %20 = zext i32 %7 to i64
  %21 = zext i32 %18 to i64
  %22 = getelementptr double, double* %4, i64 1
  %23 = getelementptr double, double* %5, i64 1
  br label %24

24:                                               ; preds = %202, %11
  %25 = phi i64 [ 0, %11 ], [ %203, %202 ]
  br i1 %12, label %26, label %202

26:                                               ; preds = %24
  %27 = getelementptr inbounds double, double* %2, i64 %25
  %28 = trunc i64 %25 to i32
  br label %30

29:                                               ; preds = %202, %9
  ret void

30:                                               ; preds = %199, %26
  %31 = phi i64 [ 0, %26 ], [ %200, %199 ]
  %32 = getelementptr inbounds double, double* %3, i64 %31
  %33 = load double, double* %32, align 8, !tbaa !4
  %34 = load double, double* %27, align 8, !tbaa !4
  %35 = fmul fast double %34, %33
  %36 = tail call fast double @llvm.cos.f64(double %35)
  %37 = fdiv fast double %36, %35
  %38 = mul nsw i64 %31, %19
  %39 = add nuw nsw i64 %38, %25
  %40 = getelementptr inbounds double, double* %0, i64 %39
  store double %37, double* %40, align 8, !tbaa !4
  %41 = load double, double* %32, align 8, !tbaa !4
  %42 = fneg fast double %36
  %43 = fmul fast double %35, %35
  %44 = fdiv fast double %42, %43
  %45 = tail call fast double @llvm.sin.f64(double %35)
  %46 = fdiv fast double %45, %35
  %47 = fsub fast double %44, %46
  %48 = fmul fast double %47, %41
  %49 = getelementptr inbounds double, double* %1, i64 %39
  store double %48, double* %49, align 8, !tbaa !4
  %50 = add nsw i64 %31, %17
  %51 = getelementptr inbounds double, double* %3, i64 %50
  %52 = load double, double* %51, align 8, !tbaa !4
  %53 = load double, double* %27, align 8, !tbaa !4
  %54 = fmul fast double %53, %52
  %55 = tail call fast double @llvm.cos.f64(double %54)
  %56 = fmul fast double %54, %54
  %57 = fdiv fast double %55, %56
  %58 = tail call fast double @llvm.sin.f64(double %54)
  %59 = fdiv fast double %58, %54
  %60 = fadd fast double %57, %59
  %61 = trunc i64 %38 to i32
  %62 = add i32 %61, %28
  %63 = add i32 %62, %13
  %64 = sext i32 %63 to i64
  %65 = getelementptr inbounds double, double* %0, i64 %64
  store double %60, double* %65, align 8, !tbaa !4
  %66 = load double, double* %51, align 8, !tbaa !4
  %67 = fdiv fast double %55, %54
  %68 = fmul fast double %56, %54
  %69 = fmul fast double %55, -2.000000e+00
  %70 = fdiv fast double %69, %68
  %71 = fadd fast double %70, %67
  %72 = fmul fast double %58, -2.000000e+00
  %73 = fdiv fast double %72, %56
  %74 = fadd fast double %71, %73
  %75 = fmul fast double %74, %66
  %76 = getelementptr inbounds double, double* %1, i64 %64
  store double %75, double* %76, align 8, !tbaa !4
  br i1 %14, label %199, label %77

77:                                               ; preds = %30, %154
  %78 = phi i64 [ %170, %154 ], [ 0, %30 ]
  %79 = phi i64 [ %155, %154 ], [ 2, %30 ]
  %80 = phi i64 [ %168, %154 ], [ 3, %30 ]
  %81 = add i64 %78, 3
  %82 = getelementptr double, double* %4, i64 %81
  %83 = getelementptr double, double* %5, i64 %81
  %84 = mul nsw i64 %79, %17
  %85 = add nsw i64 %84, %31
  %86 = getelementptr inbounds double, double* %3, i64 %85
  %87 = load double, double* %86, align 8, !tbaa !4
  %88 = load double, double* %27, align 8, !tbaa !4
  %89 = fmul fast double %88, %87
  %90 = tail call fast double @llvm.cos.f64(double %89)
  %91 = fneg fast double %90
  %92 = fdiv fast double %91, %89
  store double %92, double* %4, align 8, !tbaa !4
  %93 = fmul fast double %89, %89
  %94 = fdiv fast double %91, %93
  %95 = tail call fast double @llvm.sin.f64(double %89)
  %96 = fdiv fast double %95, %89
  %97 = fsub fast double %94, %96
  store double %97, double* %15, align 8, !tbaa !4
  %98 = load double, double* %86, align 8, !tbaa !4
  %99 = fdiv fast double %90, %93
  %100 = fadd fast double %99, %96
  %101 = fmul fast double %100, %98
  store double %101, double* %5, align 8, !tbaa !4
  %102 = load double, double* %86, align 8, !tbaa !4
  %103 = fmul fast double %90, 2.000000e+00
  %104 = fmul fast double %93, %89
  %105 = fdiv fast double %103, %104
  %106 = fdiv fast double %90, %89
  %107 = fsub fast double %105, %106
  %108 = fmul fast double %95, 2.000000e+00
  %109 = fdiv fast double %108, %93
  %110 = fadd fast double %107, %109
  %111 = fmul fast double %110, %102
  store double %111, double* %16, align 8, !tbaa !4
  %112 = icmp ugt double* %83, %4
  %113 = icmp ugt double* %82, %5
  %114 = and i1 %112, %113
  br i1 %114, label %115, label %149

115:                                              ; preds = %77
  %116 = fdiv fast double 1.000000e+00, %89
  %117 = fdiv fast double 1.000000e+00, %93
  br label %118

118:                                              ; preds = %115, %118
  %119 = phi i64 [ %120, %118 ], [ 2, %115 ]
  %120 = add nuw nsw i64 %119, 1
  %121 = trunc i64 %120 to i32
  %122 = shl nuw nsw i32 %121, 1
  %123 = add nsw i32 %122, -3
  %124 = sitofp i32 %123 to double
  %125 = fmul fast double %124, %116
  %126 = add nsw i64 %119, -1
  %127 = getelementptr inbounds double, double* %4, i64 %126
  %128 = load double, double* %127, align 8, !tbaa !4
  %129 = fmul fast double %128, %125
  %130 = add nsw i64 %119, -2
  %131 = getelementptr inbounds double, double* %4, i64 %130
  %132 = load double, double* %131, align 8, !tbaa !4
  %133 = fsub fast double %129, %132
  %134 = getelementptr inbounds double, double* %4, i64 %119
  store double %133, double* %134, align 8, !tbaa !4
  %135 = getelementptr inbounds double, double* %5, i64 %126
  %136 = load double, double* %135, align 8, !tbaa !4
  %137 = fmul fast double %136, %125
  %138 = load double, double* %86, align 8, !tbaa !4
  %139 = fneg fast double %138
  %140 = fmul fast double %139, %124
  %141 = fmul fast double %128, %140
  %142 = fmul fast double %141, %117
  %143 = getelementptr inbounds double, double* %5, i64 %130
  %144 = load double, double* %143, align 8, !tbaa !4
  %145 = fsub fast double %137, %144
  %146 = fadd fast double %145, %142
  %147 = getelementptr inbounds double, double* %5, i64 %119
  store double %146, double* %147, align 8, !tbaa !4
  %148 = icmp eq i64 %120, %80
  br i1 %148, label %154, label %118

149:                                              ; preds = %77
  %150 = load double, double* %22, align 8
  %151 = load double, double* %23, align 8
  %152 = fdiv fast double 1.000000e+00, %89
  %153 = fdiv fast double 1.000000e+00, %93
  br label %171

154:                                              ; preds = %171, %118
  %155 = add nuw nsw i64 %79, 1
  %156 = getelementptr inbounds double, double* %4, i64 %79
  %157 = load double, double* %156, align 8, !tbaa !4
  %158 = fneg fast double %157
  %159 = trunc i64 %79 to i32
  %160 = mul i32 %13, %159
  %161 = add i32 %62, %160
  %162 = sext i32 %161 to i64
  %163 = getelementptr inbounds double, double* %0, i64 %162
  store double %158, double* %163, align 8, !tbaa !4
  %164 = getelementptr inbounds double, double* %5, i64 %79
  %165 = load double, double* %164, align 8, !tbaa !4
  %166 = fneg fast double %165
  %167 = getelementptr inbounds double, double* %1, i64 %162
  store double %166, double* %167, align 8, !tbaa !4
  %168 = add nuw nsw i64 %80, 1
  %169 = icmp eq i64 %155, %21
  %170 = add i64 %78, 1
  br i1 %169, label %199, label %77

171:                                              ; preds = %171, %149
  %172 = phi double [ %151, %149 ], [ %196, %171 ]
  %173 = phi double [ %150, %149 ], [ %185, %171 ]
  %174 = phi i64 [ 2, %149 ], [ %175, %171 ]
  %175 = add nuw nsw i64 %174, 1
  %176 = trunc i64 %175 to i32
  %177 = shl nuw nsw i32 %176, 1
  %178 = add nsw i32 %177, -3
  %179 = sitofp i32 %178 to double
  %180 = fmul fast double %179, %152
  %181 = fmul fast double %173, %180
  %182 = add nsw i64 %174, -2
  %183 = getelementptr inbounds double, double* %4, i64 %182
  %184 = load double, double* %183, align 8, !tbaa !4
  %185 = fsub fast double %181, %184
  %186 = getelementptr inbounds double, double* %4, i64 %174
  store double %185, double* %186, align 8, !tbaa !4
  %187 = fmul fast double %172, %180
  %188 = load double, double* %86, align 8, !tbaa !4
  %189 = fneg fast double %188
  %190 = fmul fast double %189, %179
  %191 = fmul fast double %173, %190
  %192 = fmul fast double %191, %153
  %193 = getelementptr inbounds double, double* %5, i64 %182
  %194 = load double, double* %193, align 8, !tbaa !4
  %195 = fsub fast double %187, %194
  %196 = fadd fast double %195, %192
  %197 = getelementptr inbounds double, double* %5, i64 %174
  store double %196, double* %197, align 8, !tbaa !4
  %198 = icmp eq i64 %175, %80
  br i1 %198, label %154, label %171

199:                                              ; preds = %154, %30
  %200 = add nuw nsw i64 %31, 1
  %201 = icmp eq i64 %200, %20
  br i1 %201, label %202, label %30

202:                                              ; preds = %199, %24
  %203 = add nuw nsw i64 %25, 1
  %204 = icmp eq i64 %203, %19
  br i1 %204, label %29, label %24
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z23cpuSphericalBesselDerivIfEvPT_S1_S1_S1_S1_S1_iii(float* %0, float* %1, float* %2, float* %3, float* %4, float* %5, i32 %6, i32 %7, i32 %8) local_unnamed_addr #0 {
  %10 = icmp sgt i32 %8, 0
  br i1 %10, label %11, label %29

11:                                               ; preds = %9
  %12 = icmp sgt i32 %7, 0
  %13 = mul i32 %8, %7
  %14 = icmp slt i32 %6, 2
  %15 = getelementptr inbounds float, float* %4, i64 1
  %16 = getelementptr inbounds float, float* %5, i64 1
  %17 = sext i32 %7 to i64
  %18 = add i32 %6, 1
  %19 = zext i32 %8 to i64
  %20 = zext i32 %7 to i64
  %21 = zext i32 %18 to i64
  %22 = getelementptr float, float* %4, i64 1
  %23 = getelementptr float, float* %5, i64 1
  br label %24

24:                                               ; preds = %202, %11
  %25 = phi i64 [ 0, %11 ], [ %203, %202 ]
  br i1 %12, label %26, label %202

26:                                               ; preds = %24
  %27 = getelementptr inbounds float, float* %2, i64 %25
  %28 = trunc i64 %25 to i32
  br label %30

29:                                               ; preds = %202, %9
  ret void

30:                                               ; preds = %199, %26
  %31 = phi i64 [ 0, %26 ], [ %200, %199 ]
  %32 = getelementptr inbounds float, float* %3, i64 %31
  %33 = load float, float* %32, align 4, !tbaa !8
  %34 = load float, float* %27, align 4, !tbaa !8
  %35 = fmul fast float %34, %33
  %36 = tail call fast float @llvm.cos.f32(float %35) #8
  %37 = fdiv fast float %36, %35
  %38 = mul nsw i64 %31, %19
  %39 = add nuw nsw i64 %38, %25
  %40 = getelementptr inbounds float, float* %0, i64 %39
  store float %37, float* %40, align 4, !tbaa !8
  %41 = load float, float* %32, align 4, !tbaa !8
  %42 = fneg fast float %36
  %43 = fmul fast float %35, %35
  %44 = fdiv fast float %42, %43
  %45 = tail call fast float @llvm.sin.f32(float %35) #8
  %46 = fdiv fast float %45, %35
  %47 = fsub fast float %44, %46
  %48 = fmul fast float %47, %41
  %49 = getelementptr inbounds float, float* %1, i64 %39
  store float %48, float* %49, align 4, !tbaa !8
  %50 = add nsw i64 %31, %17
  %51 = getelementptr inbounds float, float* %3, i64 %50
  %52 = load float, float* %51, align 4, !tbaa !8
  %53 = load float, float* %27, align 4, !tbaa !8
  %54 = fmul fast float %53, %52
  %55 = tail call fast float @llvm.cos.f32(float %54) #8
  %56 = fmul fast float %54, %54
  %57 = fdiv fast float %55, %56
  %58 = tail call fast float @llvm.sin.f32(float %54) #8
  %59 = fdiv fast float %58, %54
  %60 = fadd fast float %57, %59
  %61 = trunc i64 %38 to i32
  %62 = add i32 %61, %28
  %63 = add i32 %62, %13
  %64 = sext i32 %63 to i64
  %65 = getelementptr inbounds float, float* %0, i64 %64
  store float %60, float* %65, align 4, !tbaa !8
  %66 = load float, float* %51, align 4, !tbaa !8
  %67 = fdiv fast float %55, %54
  %68 = fmul fast float %56, %54
  %69 = fmul fast float %55, -2.000000e+00
  %70 = fdiv fast float %69, %68
  %71 = fadd fast float %70, %67
  %72 = fmul fast float %58, -2.000000e+00
  %73 = fdiv fast float %72, %56
  %74 = fadd fast float %71, %73
  %75 = fmul fast float %74, %66
  %76 = getelementptr inbounds float, float* %1, i64 %64
  store float %75, float* %76, align 4, !tbaa !8
  br i1 %14, label %199, label %77

77:                                               ; preds = %30, %154
  %78 = phi i64 [ %170, %154 ], [ 0, %30 ]
  %79 = phi i64 [ %155, %154 ], [ 2, %30 ]
  %80 = phi i64 [ %168, %154 ], [ 3, %30 ]
  %81 = add i64 %78, 3
  %82 = getelementptr float, float* %4, i64 %81
  %83 = getelementptr float, float* %5, i64 %81
  %84 = mul nsw i64 %79, %17
  %85 = add nsw i64 %84, %31
  %86 = getelementptr inbounds float, float* %3, i64 %85
  %87 = load float, float* %86, align 4, !tbaa !8
  %88 = load float, float* %27, align 4, !tbaa !8
  %89 = fmul fast float %88, %87
  %90 = tail call fast float @llvm.cos.f32(float %89) #8
  %91 = fneg fast float %90
  %92 = fdiv fast float %91, %89
  store float %92, float* %4, align 4, !tbaa !8
  %93 = fmul fast float %89, %89
  %94 = fdiv fast float %91, %93
  %95 = tail call fast float @llvm.sin.f32(float %89) #8
  %96 = fdiv fast float %95, %89
  %97 = fsub fast float %94, %96
  store float %97, float* %15, align 4, !tbaa !8
  %98 = load float, float* %86, align 4, !tbaa !8
  %99 = fdiv fast float %90, %93
  %100 = fadd fast float %99, %96
  %101 = fmul fast float %100, %98
  store float %101, float* %5, align 4, !tbaa !8
  %102 = load float, float* %86, align 4, !tbaa !8
  %103 = fmul fast float %90, 2.000000e+00
  %104 = fmul fast float %93, %89
  %105 = fdiv fast float %103, %104
  %106 = fdiv fast float %90, %89
  %107 = fsub fast float %105, %106
  %108 = fmul fast float %95, 2.000000e+00
  %109 = fdiv fast float %108, %93
  %110 = fadd fast float %107, %109
  %111 = fmul fast float %110, %102
  store float %111, float* %16, align 4, !tbaa !8
  %112 = icmp ugt float* %83, %4
  %113 = icmp ugt float* %82, %5
  %114 = and i1 %112, %113
  br i1 %114, label %115, label %149

115:                                              ; preds = %77
  %116 = fdiv fast float 1.000000e+00, %89
  %117 = fdiv fast float 1.000000e+00, %93
  br label %118

118:                                              ; preds = %115, %118
  %119 = phi i64 [ %120, %118 ], [ 2, %115 ]
  %120 = add nuw nsw i64 %119, 1
  %121 = trunc i64 %120 to i32
  %122 = shl nuw nsw i32 %121, 1
  %123 = add nsw i32 %122, -3
  %124 = sitofp i32 %123 to float
  %125 = fmul fast float %124, %116
  %126 = add nsw i64 %119, -1
  %127 = getelementptr inbounds float, float* %4, i64 %126
  %128 = load float, float* %127, align 4, !tbaa !8
  %129 = fmul fast float %128, %125
  %130 = add nsw i64 %119, -2
  %131 = getelementptr inbounds float, float* %4, i64 %130
  %132 = load float, float* %131, align 4, !tbaa !8
  %133 = fsub fast float %129, %132
  %134 = getelementptr inbounds float, float* %4, i64 %119
  store float %133, float* %134, align 4, !tbaa !8
  %135 = getelementptr inbounds float, float* %5, i64 %126
  %136 = load float, float* %135, align 4, !tbaa !8
  %137 = fmul fast float %136, %125
  %138 = load float, float* %86, align 4, !tbaa !8
  %139 = fneg fast float %138
  %140 = fmul fast float %139, %124
  %141 = fmul fast float %128, %140
  %142 = fmul fast float %141, %117
  %143 = getelementptr inbounds float, float* %5, i64 %130
  %144 = load float, float* %143, align 4, !tbaa !8
  %145 = fsub fast float %137, %144
  %146 = fadd fast float %145, %142
  %147 = getelementptr inbounds float, float* %5, i64 %119
  store float %146, float* %147, align 4, !tbaa !8
  %148 = icmp eq i64 %120, %80
  br i1 %148, label %154, label %118

149:                                              ; preds = %77
  %150 = load float, float* %22, align 4
  %151 = load float, float* %23, align 4
  %152 = fdiv fast float 1.000000e+00, %89
  %153 = fdiv fast float 1.000000e+00, %93
  br label %171

154:                                              ; preds = %171, %118
  %155 = add nuw nsw i64 %79, 1
  %156 = getelementptr inbounds float, float* %4, i64 %79
  %157 = load float, float* %156, align 4, !tbaa !8
  %158 = fneg fast float %157
  %159 = trunc i64 %79 to i32
  %160 = mul i32 %13, %159
  %161 = add i32 %62, %160
  %162 = sext i32 %161 to i64
  %163 = getelementptr inbounds float, float* %0, i64 %162
  store float %158, float* %163, align 4, !tbaa !8
  %164 = getelementptr inbounds float, float* %5, i64 %79
  %165 = load float, float* %164, align 4, !tbaa !8
  %166 = fneg fast float %165
  %167 = getelementptr inbounds float, float* %1, i64 %162
  store float %166, float* %167, align 4, !tbaa !8
  %168 = add nuw nsw i64 %80, 1
  %169 = icmp eq i64 %155, %21
  %170 = add i64 %78, 1
  br i1 %169, label %199, label %77

171:                                              ; preds = %171, %149
  %172 = phi float [ %151, %149 ], [ %196, %171 ]
  %173 = phi float [ %150, %149 ], [ %185, %171 ]
  %174 = phi i64 [ 2, %149 ], [ %175, %171 ]
  %175 = add nuw nsw i64 %174, 1
  %176 = trunc i64 %175 to i32
  %177 = shl nuw nsw i32 %176, 1
  %178 = add nsw i32 %177, -3
  %179 = sitofp i32 %178 to float
  %180 = fmul fast float %179, %152
  %181 = fmul fast float %173, %180
  %182 = add nsw i64 %174, -2
  %183 = getelementptr inbounds float, float* %4, i64 %182
  %184 = load float, float* %183, align 4, !tbaa !8
  %185 = fsub fast float %181, %184
  %186 = getelementptr inbounds float, float* %4, i64 %174
  store float %185, float* %186, align 4, !tbaa !8
  %187 = fmul fast float %172, %180
  %188 = load float, float* %86, align 4, !tbaa !8
  %189 = fneg fast float %188
  %190 = fmul fast float %189, %179
  %191 = fmul fast float %173, %190
  %192 = fmul fast float %191, %153
  %193 = getelementptr inbounds float, float* %5, i64 %182
  %194 = load float, float* %193, align 4, !tbaa !8
  %195 = fsub fast float %187, %194
  %196 = fadd fast float %195, %192
  %197 = getelementptr inbounds float, float* %5, i64 %174
  store float %196, float* %197, align 4, !tbaa !8
  %198 = icmp eq i64 %175, %80
  br i1 %198, label %154, label %171

199:                                              ; preds = %154, %30
  %200 = add nuw nsw i64 %31, 1
  %201 = icmp eq i64 %200, %20
  br i1 %201, label %202, label %30

202:                                              ; preds = %199, %24
  %203 = add nuw nsw i64 %25, 1
  %204 = icmp eq i64 %203, %19
  br i1 %204, label %29, label %24
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z27cpuRadialSphericalHarmonicsIdEvPT_S1_S1_S1_S1_iii(double* %0, double* %1, double* %2, double* %3, double* %4, i32 %5, i32 %6, i32 %7) local_unnamed_addr #0 {
  %9 = icmp sgt i32 %7, 0
  br i1 %9, label %10, label %70

10:                                               ; preds = %8
  %11 = icmp sgt i32 %6, 0
  %12 = icmp slt i32 %5, 0
  %13 = mul i32 %7, %6
  br i1 %11, label %14, label %70

14:                                               ; preds = %10
  %15 = zext i32 %7 to i64
  %16 = add i32 %5, 1
  %17 = zext i32 %16 to i64
  br label %18

18:                                               ; preds = %67, %14
  %19 = phi i64 [ 0, %14 ], [ %68, %67 ]
  br i1 %12, label %67, label %20

20:                                               ; preds = %18
  %21 = trunc i64 %19 to i32
  br label %60

22:                                               ; preds = %22, %43
  %23 = phi i64 [ %41, %22 ], [ 0, %43 ]
  %24 = load double, double* %53, align 8, !tbaa !4
  %25 = add nuw nsw i64 %23, %56
  %26 = mul nsw i64 %25, %15
  %27 = add nuw nsw i64 %26, %19
  %28 = getelementptr inbounds double, double* %3, i64 %27
  %29 = load double, double* %28, align 8, !tbaa !4
  %30 = fmul fast double %29, %24
  %31 = trunc i64 %26 to i32
  %32 = mul i32 %31, %6
  %33 = add i32 %63, %32
  %34 = sext i32 %33 to i64
  %35 = getelementptr inbounds double, double* %0, i64 %34
  store double %30, double* %35, align 8, !tbaa !4
  %36 = load double, double* %53, align 8, !tbaa !4
  %37 = getelementptr inbounds double, double* %4, i64 %27
  %38 = load double, double* %37, align 8, !tbaa !4
  %39 = fmul fast double %38, %36
  %40 = getelementptr inbounds double, double* %1, i64 %34
  store double %39, double* %40, align 8, !tbaa !4
  %41 = add nuw nsw i64 %23, 1
  %42 = icmp eq i64 %41, %45
  br i1 %42, label %57, label %22

43:                                               ; preds = %57, %60
  %44 = phi i64 [ %47, %57 ], [ 0, %60 ]
  %45 = phi i64 [ %58, %57 ], [ 1, %60 ]
  %46 = phi i32 [ %48, %57 ], [ 0, %60 ]
  %47 = add nuw nsw i64 %44, 1
  %48 = add nuw nsw i32 %46, 1
  %49 = trunc i64 %44 to i32
  %50 = mul i32 %13, %49
  %51 = add i32 %63, %50
  %52 = sext i32 %51 to i64
  %53 = getelementptr inbounds double, double* %2, i64 %52
  %54 = mul nsw i32 %48, %49
  %55 = lshr i32 %54, 1
  %56 = zext i32 %55 to i64
  br label %22

57:                                               ; preds = %22
  %58 = add nuw nsw i64 %45, 1
  %59 = icmp eq i64 %47, %17
  br i1 %59, label %64, label %43

60:                                               ; preds = %64, %20
  %61 = phi i32 [ %65, %64 ], [ 0, %20 ]
  %62 = mul nsw i32 %61, %7
  %63 = add i32 %62, %21
  br label %43

64:                                               ; preds = %57
  %65 = add nuw nsw i32 %61, 1
  %66 = icmp eq i32 %65, %6
  br i1 %66, label %67, label %60

67:                                               ; preds = %64, %18
  %68 = add nuw nsw i64 %19, 1
  %69 = icmp eq i64 %68, %15
  br i1 %69, label %70, label %18

70:                                               ; preds = %67, %10, %8
  ret void
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z27cpuRadialSphericalHarmonicsIfEvPT_S1_S1_S1_S1_iii(float* %0, float* %1, float* %2, float* %3, float* %4, i32 %5, i32 %6, i32 %7) local_unnamed_addr #0 {
  %9 = icmp sgt i32 %7, 0
  br i1 %9, label %10, label %70

10:                                               ; preds = %8
  %11 = icmp sgt i32 %6, 0
  %12 = icmp slt i32 %5, 0
  %13 = mul i32 %7, %6
  br i1 %11, label %14, label %70

14:                                               ; preds = %10
  %15 = zext i32 %7 to i64
  %16 = add i32 %5, 1
  %17 = zext i32 %16 to i64
  br label %18

18:                                               ; preds = %67, %14
  %19 = phi i64 [ 0, %14 ], [ %68, %67 ]
  br i1 %12, label %67, label %20

20:                                               ; preds = %18
  %21 = trunc i64 %19 to i32
  br label %60

22:                                               ; preds = %22, %43
  %23 = phi i64 [ %41, %22 ], [ 0, %43 ]
  %24 = load float, float* %53, align 4, !tbaa !8
  %25 = add nuw nsw i64 %23, %56
  %26 = mul nsw i64 %25, %15
  %27 = add nuw nsw i64 %26, %19
  %28 = getelementptr inbounds float, float* %3, i64 %27
  %29 = load float, float* %28, align 4, !tbaa !8
  %30 = fmul fast float %29, %24
  %31 = trunc i64 %26 to i32
  %32 = mul i32 %31, %6
  %33 = add i32 %63, %32
  %34 = sext i32 %33 to i64
  %35 = getelementptr inbounds float, float* %0, i64 %34
  store float %30, float* %35, align 4, !tbaa !8
  %36 = load float, float* %53, align 4, !tbaa !8
  %37 = getelementptr inbounds float, float* %4, i64 %27
  %38 = load float, float* %37, align 4, !tbaa !8
  %39 = fmul fast float %38, %36
  %40 = getelementptr inbounds float, float* %1, i64 %34
  store float %39, float* %40, align 4, !tbaa !8
  %41 = add nuw nsw i64 %23, 1
  %42 = icmp eq i64 %41, %45
  br i1 %42, label %57, label %22

43:                                               ; preds = %57, %60
  %44 = phi i64 [ %47, %57 ], [ 0, %60 ]
  %45 = phi i64 [ %58, %57 ], [ 1, %60 ]
  %46 = phi i32 [ %48, %57 ], [ 0, %60 ]
  %47 = add nuw nsw i64 %44, 1
  %48 = add nuw nsw i32 %46, 1
  %49 = trunc i64 %44 to i32
  %50 = mul i32 %13, %49
  %51 = add i32 %63, %50
  %52 = sext i32 %51 to i64
  %53 = getelementptr inbounds float, float* %2, i64 %52
  %54 = mul nsw i32 %48, %49
  %55 = lshr i32 %54, 1
  %56 = zext i32 %55 to i64
  br label %22

57:                                               ; preds = %22
  %58 = add nuw nsw i64 %45, 1
  %59 = icmp eq i64 %47, %17
  br i1 %59, label %64, label %43

60:                                               ; preds = %64, %20
  %61 = phi i32 [ %65, %64 ], [ 0, %20 ]
  %62 = mul nsw i32 %61, %7
  %63 = add i32 %62, %21
  br label %43

64:                                               ; preds = %57
  %65 = add nuw nsw i32 %61, 1
  %66 = icmp eq i32 %65, %6
  br i1 %66, label %67, label %60

67:                                               ; preds = %64, %18
  %68 = add nuw nsw i64 %19, 1
  %69 = icmp eq i64 %68, %15
  br i1 %69, label %70, label %18

70:                                               ; preds = %67, %10, %8
  ret void
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z35cpuWeightedRadialSphericalHarmonicsIdEvPT_S1_S1_S1_S1_S1_iii(double* %0, double* %1, double* %2, double* %3, double* %4, double* %5, i32 %6, i32 %7, i32 %8) local_unnamed_addr #0 {
  %10 = icmp sgt i32 %8, 0
  br i1 %10, label %11, label %76

11:                                               ; preds = %9
  %12 = icmp sgt i32 %7, 0
  %13 = icmp slt i32 %6, 0
  %14 = mul i32 %8, %7
  br i1 %12, label %15, label %76

15:                                               ; preds = %11
  %16 = zext i32 %8 to i64
  %17 = add i32 %6, 1
  %18 = zext i32 %17 to i64
  br label %19

19:                                               ; preds = %73, %15
  %20 = phi i64 [ 0, %15 ], [ %74, %73 ]
  %21 = getelementptr inbounds double, double* %5, i64 %20
  br i1 %13, label %73, label %22

22:                                               ; preds = %19
  %23 = trunc i64 %20 to i32
  br label %66

24:                                               ; preds = %24, %49
  %25 = phi i64 [ %47, %24 ], [ 0, %49 ]
  %26 = load double, double* %21, align 8, !tbaa !4
  %27 = load double, double* %59, align 8, !tbaa !4
  %28 = fmul fast double %27, %26
  %29 = add nuw nsw i64 %25, %62
  %30 = mul nsw i64 %29, %16
  %31 = add nuw nsw i64 %30, %20
  %32 = getelementptr inbounds double, double* %3, i64 %31
  %33 = load double, double* %32, align 8, !tbaa !4
  %34 = fmul fast double %28, %33
  %35 = trunc i64 %30 to i32
  %36 = mul i32 %35, %7
  %37 = add i32 %69, %36
  %38 = sext i32 %37 to i64
  %39 = getelementptr inbounds double, double* %0, i64 %38
  store double %34, double* %39, align 8, !tbaa !4
  %40 = load double, double* %21, align 8, !tbaa !4
  %41 = load double, double* %59, align 8, !tbaa !4
  %42 = fmul fast double %41, %40
  %43 = getelementptr inbounds double, double* %4, i64 %31
  %44 = load double, double* %43, align 8, !tbaa !4
  %45 = fmul fast double %42, %44
  %46 = getelementptr inbounds double, double* %1, i64 %38
  store double %45, double* %46, align 8, !tbaa !4
  %47 = add nuw nsw i64 %25, 1
  %48 = icmp eq i64 %47, %51
  br i1 %48, label %63, label %24

49:                                               ; preds = %63, %66
  %50 = phi i64 [ %53, %63 ], [ 0, %66 ]
  %51 = phi i64 [ %64, %63 ], [ 1, %66 ]
  %52 = phi i32 [ %54, %63 ], [ 0, %66 ]
  %53 = add nuw nsw i64 %50, 1
  %54 = add nuw nsw i32 %52, 1
  %55 = trunc i64 %50 to i32
  %56 = mul i32 %14, %55
  %57 = add i32 %69, %56
  %58 = sext i32 %57 to i64
  %59 = getelementptr inbounds double, double* %2, i64 %58
  %60 = mul nsw i32 %54, %55
  %61 = lshr i32 %60, 1
  %62 = zext i32 %61 to i64
  br label %24

63:                                               ; preds = %24
  %64 = add nuw nsw i64 %51, 1
  %65 = icmp eq i64 %53, %18
  br i1 %65, label %70, label %49

66:                                               ; preds = %70, %22
  %67 = phi i32 [ %71, %70 ], [ 0, %22 ]
  %68 = mul nsw i32 %67, %8
  %69 = add i32 %68, %23
  br label %49

70:                                               ; preds = %63
  %71 = add nuw nsw i32 %67, 1
  %72 = icmp eq i32 %71, %7
  br i1 %72, label %73, label %66

73:                                               ; preds = %70, %19
  %74 = add nuw nsw i64 %20, 1
  %75 = icmp eq i64 %74, %16
  br i1 %75, label %76, label %19

76:                                               ; preds = %73, %11, %9
  ret void
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z35cpuWeightedRadialSphericalHarmonicsIfEvPT_S1_S1_S1_S1_S1_iii(float* %0, float* %1, float* %2, float* %3, float* %4, float* %5, i32 %6, i32 %7, i32 %8) local_unnamed_addr #0 {
  %10 = icmp sgt i32 %8, 0
  br i1 %10, label %11, label %76

11:                                               ; preds = %9
  %12 = icmp sgt i32 %7, 0
  %13 = icmp slt i32 %6, 0
  %14 = mul i32 %8, %7
  br i1 %12, label %15, label %76

15:                                               ; preds = %11
  %16 = zext i32 %8 to i64
  %17 = add i32 %6, 1
  %18 = zext i32 %17 to i64
  br label %19

19:                                               ; preds = %73, %15
  %20 = phi i64 [ 0, %15 ], [ %74, %73 ]
  %21 = getelementptr inbounds float, float* %5, i64 %20
  br i1 %13, label %73, label %22

22:                                               ; preds = %19
  %23 = trunc i64 %20 to i32
  br label %66

24:                                               ; preds = %24, %49
  %25 = phi i64 [ %47, %24 ], [ 0, %49 ]
  %26 = load float, float* %21, align 4, !tbaa !8
  %27 = load float, float* %59, align 4, !tbaa !8
  %28 = fmul fast float %27, %26
  %29 = add nuw nsw i64 %25, %62
  %30 = mul nsw i64 %29, %16
  %31 = add nuw nsw i64 %30, %20
  %32 = getelementptr inbounds float, float* %3, i64 %31
  %33 = load float, float* %32, align 4, !tbaa !8
  %34 = fmul fast float %28, %33
  %35 = trunc i64 %30 to i32
  %36 = mul i32 %35, %7
  %37 = add i32 %69, %36
  %38 = sext i32 %37 to i64
  %39 = getelementptr inbounds float, float* %0, i64 %38
  store float %34, float* %39, align 4, !tbaa !8
  %40 = load float, float* %21, align 4, !tbaa !8
  %41 = load float, float* %59, align 4, !tbaa !8
  %42 = fmul fast float %41, %40
  %43 = getelementptr inbounds float, float* %4, i64 %31
  %44 = load float, float* %43, align 4, !tbaa !8
  %45 = fmul fast float %42, %44
  %46 = getelementptr inbounds float, float* %1, i64 %38
  store float %45, float* %46, align 4, !tbaa !8
  %47 = add nuw nsw i64 %25, 1
  %48 = icmp eq i64 %47, %51
  br i1 %48, label %63, label %24

49:                                               ; preds = %63, %66
  %50 = phi i64 [ %53, %63 ], [ 0, %66 ]
  %51 = phi i64 [ %64, %63 ], [ 1, %66 ]
  %52 = phi i32 [ %54, %63 ], [ 0, %66 ]
  %53 = add nuw nsw i64 %50, 1
  %54 = add nuw nsw i32 %52, 1
  %55 = trunc i64 %50 to i32
  %56 = mul i32 %14, %55
  %57 = add i32 %69, %56
  %58 = sext i32 %57 to i64
  %59 = getelementptr inbounds float, float* %2, i64 %58
  %60 = mul nsw i32 %54, %55
  %61 = lshr i32 %60, 1
  %62 = zext i32 %61 to i64
  br label %24

63:                                               ; preds = %24
  %64 = add nuw nsw i64 %51, 1
  %65 = icmp eq i64 %53, %18
  br i1 %65, label %70, label %49

66:                                               ; preds = %70, %22
  %67 = phi i32 [ %71, %70 ], [ 0, %22 ]
  %68 = mul nsw i32 %67, %8
  %69 = add i32 %68, %23
  br label %49

70:                                               ; preds = %63
  %71 = add nuw nsw i32 %67, 1
  %72 = icmp eq i32 %71, %7
  br i1 %72, label %73, label %66

73:                                               ; preds = %70, %19
  %74 = add nuw nsw i64 %20, 1
  %75 = icmp eq i64 %74, %16
  br i1 %75, label %76, label %19

76:                                               ; preds = %73, %11, %9
  ret void
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z30cpuRadialSphericalHarmonicsSumIdEvPT_S1_S1_S1_S1_iii(double* %0, double* %1, double* %2, double* %3, double* %4, i32 %5, i32 %6, i32 %7) local_unnamed_addr #0 {
  %9 = icmp sgt i32 %6, 0
  br i1 %9, label %10, label %67

10:                                               ; preds = %8
  %11 = icmp slt i32 %5, 0
  %12 = icmp sgt i32 %7, 0
  %13 = zext i32 %6 to i64
  %14 = sext i32 %7 to i64
  %15 = zext i32 %7 to i64
  %16 = icmp eq i32 %6, 1
  br label %17

17:                                               ; preds = %159, %10
  %18 = phi i64 [ 0, %10 ], [ %160, %159 ]
  br i1 %11, label %159, label %19

19:                                               ; preds = %17
  br i1 %12, label %20, label %72

20:                                               ; preds = %19
  %21 = trunc i64 %18 to i32
  br label %22

22:                                               ; preds = %64, %20
  %23 = phi i64 [ 1, %20 ], [ %65, %64 ]
  %24 = phi i32 [ 0, %20 ], [ %25, %64 ]
  %25 = add nuw i32 %24, 1
  %26 = mul nsw i32 %25, %24
  %27 = lshr i32 %26, 1
  %28 = mul i32 %24, %6
  %29 = add i32 %28, %21
  %30 = mul i32 %29, %7
  %31 = zext i32 %27 to i64
  br label %32

32:                                               ; preds = %40, %22
  %33 = phi i64 [ %41, %40 ], [ 0, %22 ]
  %34 = add nuw nsw i64 %33, %31
  %35 = mul nsw i64 %34, %13
  %36 = add nuw nsw i64 %35, %18
  %37 = getelementptr inbounds double, double* %0, i64 %36
  store double 0.000000e+00, double* %37, align 8, !tbaa !4
  %38 = getelementptr inbounds double, double* %1, i64 %36
  store double 0.000000e+00, double* %38, align 8, !tbaa !4
  %39 = mul nsw i64 %34, %14
  br label %43

40:                                               ; preds = %43
  %41 = add nuw nsw i64 %33, 1
  %42 = icmp eq i64 %41, %23
  br i1 %42, label %64, label %32

43:                                               ; preds = %43, %32
  %44 = phi i64 [ %62, %43 ], [ 0, %32 ]
  %45 = trunc i64 %44 to i32
  %46 = add i32 %30, %45
  %47 = sext i32 %46 to i64
  %48 = getelementptr inbounds double, double* %2, i64 %47
  %49 = load double, double* %48, align 8, !tbaa !4
  %50 = add nsw i64 %44, %39
  %51 = getelementptr inbounds double, double* %3, i64 %50
  %52 = load double, double* %51, align 8, !tbaa !4
  %53 = fmul fast double %52, %49
  %54 = load double, double* %37, align 8, !tbaa !4
  %55 = fadd fast double %54, %53
  store double %55, double* %37, align 8, !tbaa !4
  %56 = load double, double* %48, align 8, !tbaa !4
  %57 = getelementptr inbounds double, double* %4, i64 %50
  %58 = load double, double* %57, align 8, !tbaa !4
  %59 = fmul fast double %58, %56
  %60 = load double, double* %38, align 8, !tbaa !4
  %61 = fadd fast double %60, %59
  store double %61, double* %38, align 8, !tbaa !4
  %62 = add nuw nsw i64 %44, 1
  %63 = icmp eq i64 %62, %15
  br i1 %63, label %40, label %43

64:                                               ; preds = %40
  %65 = add nuw nsw i64 %23, 1
  %66 = icmp eq i32 %24, %5
  br i1 %66, label %159, label %22

67:                                               ; preds = %159, %8
  ret void

68:                                               ; preds = %104, %162, %157
  %69 = add nuw nsw i64 %74, 1
  %70 = icmp eq i32 %75, %5
  %71 = add i64 %73, 1
  br i1 %70, label %159, label %72

72:                                               ; preds = %19, %68
  %73 = phi i64 [ %71, %68 ], [ 0, %19 ]
  %74 = phi i64 [ %69, %68 ], [ 1, %19 ]
  %75 = phi i32 [ %81, %68 ], [ 0, %19 ]
  %76 = add i64 %73, -3
  %77 = lshr i64 %76, 2
  %78 = add nuw nsw i64 %77, 1
  %79 = add i64 %73, 1
  %80 = add i64 %73, 1
  %81 = add nuw i32 %75, 1
  %82 = mul nsw i32 %81, %75
  %83 = lshr i32 %82, 1
  %84 = zext i32 %83 to i64
  %85 = icmp ugt i64 %80, 3
  %86 = and i1 %85, %16
  br i1 %86, label %87, label %93

87:                                               ; preds = %72
  %88 = getelementptr double, double* %1, i64 %79
  %89 = getelementptr double, double* %0, i64 %79
  %90 = icmp ugt double* %88, %0
  %91 = icmp ugt double* %89, %1
  %92 = and i1 %90, %91
  br i1 %92, label %93, label %107

93:                                               ; preds = %72, %157, %87
  %94 = phi i64 [ 0, %87 ], [ 0, %72 ], [ %108, %157 ]
  %95 = and i64 %73, 1
  %96 = icmp eq i64 %95, 0
  br i1 %96, label %97, label %104

97:                                               ; preds = %93
  %98 = add nuw nsw i64 %94, %84
  %99 = mul nsw i64 %98, %13
  %100 = add nuw nsw i64 %99, %18
  %101 = getelementptr inbounds double, double* %0, i64 %100
  store double 0.000000e+00, double* %101, align 8, !tbaa !4
  %102 = getelementptr inbounds double, double* %1, i64 %100
  store double 0.000000e+00, double* %102, align 8, !tbaa !4
  %103 = or i64 %94, 1
  br label %104

104:                                              ; preds = %97, %93
  %105 = phi i64 [ %103, %97 ], [ %94, %93 ]
  %106 = icmp eq i64 %73, %94
  br i1 %106, label %68, label %162

107:                                              ; preds = %87
  %108 = and i64 %80, -4
  %109 = and i64 %78, 1
  %110 = icmp eq i64 %77, 0
  br i1 %110, label %142, label %111

111:                                              ; preds = %107
  %112 = sub nuw nsw i64 %78, %109
  br label %113

113:                                              ; preds = %113, %111
  %114 = phi i64 [ 0, %111 ], [ %139, %113 ]
  %115 = phi i64 [ %112, %111 ], [ %140, %113 ]
  %116 = add nuw nsw i64 %114, %84
  %117 = mul nsw i64 %116, %13
  %118 = add nuw nsw i64 %117, %18
  %119 = getelementptr inbounds double, double* %0, i64 %118
  %120 = bitcast double* %119 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %120, align 8, !tbaa !4, !alias.scope !44, !noalias !47
  %121 = getelementptr inbounds double, double* %119, i64 2
  %122 = bitcast double* %121 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %122, align 8, !tbaa !4, !alias.scope !44, !noalias !47
  %123 = getelementptr inbounds double, double* %1, i64 %118
  %124 = bitcast double* %123 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %124, align 8, !tbaa !4, !alias.scope !47
  %125 = getelementptr inbounds double, double* %123, i64 2
  %126 = bitcast double* %125 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %126, align 8, !tbaa !4, !alias.scope !47
  %127 = or i64 %114, 4
  %128 = add nuw nsw i64 %127, %84
  %129 = mul nsw i64 %128, %13
  %130 = add nuw nsw i64 %129, %18
  %131 = getelementptr inbounds double, double* %0, i64 %130
  %132 = bitcast double* %131 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %132, align 8, !tbaa !4, !alias.scope !44, !noalias !47
  %133 = getelementptr inbounds double, double* %131, i64 2
  %134 = bitcast double* %133 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %134, align 8, !tbaa !4, !alias.scope !44, !noalias !47
  %135 = getelementptr inbounds double, double* %1, i64 %130
  %136 = bitcast double* %135 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %136, align 8, !tbaa !4, !alias.scope !47
  %137 = getelementptr inbounds double, double* %135, i64 2
  %138 = bitcast double* %137 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %138, align 8, !tbaa !4, !alias.scope !47
  %139 = add i64 %114, 8
  %140 = add i64 %115, -2
  %141 = icmp eq i64 %140, 0
  br i1 %141, label %142, label %113, !llvm.loop !49

142:                                              ; preds = %113, %107
  %143 = phi i64 [ 0, %107 ], [ %139, %113 ]
  %144 = icmp eq i64 %109, 0
  br i1 %144, label %157, label %145

145:                                              ; preds = %142
  %146 = add nuw nsw i64 %143, %84
  %147 = mul nsw i64 %146, %13
  %148 = add nuw nsw i64 %147, %18
  %149 = getelementptr inbounds double, double* %0, i64 %148
  %150 = bitcast double* %149 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %150, align 8, !tbaa !4, !alias.scope !44, !noalias !47
  %151 = getelementptr inbounds double, double* %149, i64 2
  %152 = bitcast double* %151 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %152, align 8, !tbaa !4, !alias.scope !44, !noalias !47
  %153 = getelementptr inbounds double, double* %1, i64 %148
  %154 = bitcast double* %153 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %154, align 8, !tbaa !4, !alias.scope !47
  %155 = getelementptr inbounds double, double* %153, i64 2
  %156 = bitcast double* %155 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %156, align 8, !tbaa !4, !alias.scope !47
  br label %157

157:                                              ; preds = %142, %145
  %158 = icmp eq i64 %80, %108
  br i1 %158, label %68, label %93

159:                                              ; preds = %68, %64, %17
  %160 = add nuw nsw i64 %18, 1
  %161 = icmp eq i64 %160, %13
  br i1 %161, label %67, label %17

162:                                              ; preds = %104, %162
  %163 = phi i64 [ %175, %162 ], [ %105, %104 ]
  %164 = add nuw nsw i64 %163, %84
  %165 = mul nsw i64 %164, %13
  %166 = add nuw nsw i64 %165, %18
  %167 = getelementptr inbounds double, double* %0, i64 %166
  store double 0.000000e+00, double* %167, align 8, !tbaa !4
  %168 = getelementptr inbounds double, double* %1, i64 %166
  store double 0.000000e+00, double* %168, align 8, !tbaa !4
  %169 = add nuw nsw i64 %163, 1
  %170 = add nuw nsw i64 %169, %84
  %171 = mul nsw i64 %170, %13
  %172 = add nuw nsw i64 %171, %18
  %173 = getelementptr inbounds double, double* %0, i64 %172
  store double 0.000000e+00, double* %173, align 8, !tbaa !4
  %174 = getelementptr inbounds double, double* %1, i64 %172
  store double 0.000000e+00, double* %174, align 8, !tbaa !4
  %175 = add nuw nsw i64 %163, 2
  %176 = icmp eq i64 %175, %74
  br i1 %176, label %68, label %162, !llvm.loop !50
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z30cpuRadialSphericalHarmonicsSumIfEvPT_S1_S1_S1_S1_iii(float* %0, float* %1, float* %2, float* %3, float* %4, i32 %5, i32 %6, i32 %7) local_unnamed_addr #0 {
  %9 = icmp sgt i32 %6, 0
  br i1 %9, label %10, label %67

10:                                               ; preds = %8
  %11 = icmp slt i32 %5, 0
  %12 = icmp sgt i32 %7, 0
  %13 = zext i32 %6 to i64
  %14 = sext i32 %7 to i64
  %15 = zext i32 %7 to i64
  %16 = icmp eq i32 %6, 1
  br label %17

17:                                               ; preds = %159, %10
  %18 = phi i64 [ 0, %10 ], [ %160, %159 ]
  br i1 %11, label %159, label %19

19:                                               ; preds = %17
  br i1 %12, label %20, label %72

20:                                               ; preds = %19
  %21 = trunc i64 %18 to i32
  br label %22

22:                                               ; preds = %64, %20
  %23 = phi i64 [ 1, %20 ], [ %65, %64 ]
  %24 = phi i32 [ 0, %20 ], [ %25, %64 ]
  %25 = add nuw i32 %24, 1
  %26 = mul nsw i32 %25, %24
  %27 = lshr i32 %26, 1
  %28 = mul i32 %24, %6
  %29 = add i32 %28, %21
  %30 = mul i32 %29, %7
  %31 = zext i32 %27 to i64
  br label %32

32:                                               ; preds = %40, %22
  %33 = phi i64 [ %41, %40 ], [ 0, %22 ]
  %34 = add nuw nsw i64 %33, %31
  %35 = mul nsw i64 %34, %13
  %36 = add nuw nsw i64 %35, %18
  %37 = getelementptr inbounds float, float* %0, i64 %36
  store float 0.000000e+00, float* %37, align 4, !tbaa !8
  %38 = getelementptr inbounds float, float* %1, i64 %36
  store float 0.000000e+00, float* %38, align 4, !tbaa !8
  %39 = mul nsw i64 %34, %14
  br label %43

40:                                               ; preds = %43
  %41 = add nuw nsw i64 %33, 1
  %42 = icmp eq i64 %41, %23
  br i1 %42, label %64, label %32

43:                                               ; preds = %43, %32
  %44 = phi i64 [ %62, %43 ], [ 0, %32 ]
  %45 = trunc i64 %44 to i32
  %46 = add i32 %30, %45
  %47 = sext i32 %46 to i64
  %48 = getelementptr inbounds float, float* %2, i64 %47
  %49 = load float, float* %48, align 4, !tbaa !8
  %50 = add nsw i64 %44, %39
  %51 = getelementptr inbounds float, float* %3, i64 %50
  %52 = load float, float* %51, align 4, !tbaa !8
  %53 = fmul fast float %52, %49
  %54 = load float, float* %37, align 4, !tbaa !8
  %55 = fadd fast float %54, %53
  store float %55, float* %37, align 4, !tbaa !8
  %56 = load float, float* %48, align 4, !tbaa !8
  %57 = getelementptr inbounds float, float* %4, i64 %50
  %58 = load float, float* %57, align 4, !tbaa !8
  %59 = fmul fast float %58, %56
  %60 = load float, float* %38, align 4, !tbaa !8
  %61 = fadd fast float %60, %59
  store float %61, float* %38, align 4, !tbaa !8
  %62 = add nuw nsw i64 %44, 1
  %63 = icmp eq i64 %62, %15
  br i1 %63, label %40, label %43

64:                                               ; preds = %40
  %65 = add nuw nsw i64 %23, 1
  %66 = icmp eq i32 %24, %5
  br i1 %66, label %159, label %22

67:                                               ; preds = %159, %8
  ret void

68:                                               ; preds = %104, %162, %157
  %69 = add nuw nsw i64 %74, 1
  %70 = icmp eq i32 %75, %5
  %71 = add i64 %73, 1
  br i1 %70, label %159, label %72

72:                                               ; preds = %19, %68
  %73 = phi i64 [ %71, %68 ], [ 0, %19 ]
  %74 = phi i64 [ %69, %68 ], [ 1, %19 ]
  %75 = phi i32 [ %81, %68 ], [ 0, %19 ]
  %76 = add i64 %73, -7
  %77 = lshr i64 %76, 3
  %78 = add nuw nsw i64 %77, 1
  %79 = add i64 %73, 1
  %80 = add i64 %73, 1
  %81 = add nuw i32 %75, 1
  %82 = mul nsw i32 %81, %75
  %83 = lshr i32 %82, 1
  %84 = zext i32 %83 to i64
  %85 = icmp ugt i64 %80, 7
  %86 = and i1 %85, %16
  br i1 %86, label %87, label %93

87:                                               ; preds = %72
  %88 = getelementptr float, float* %1, i64 %79
  %89 = getelementptr float, float* %0, i64 %79
  %90 = icmp ugt float* %88, %0
  %91 = icmp ugt float* %89, %1
  %92 = and i1 %90, %91
  br i1 %92, label %93, label %107

93:                                               ; preds = %72, %157, %87
  %94 = phi i64 [ 0, %87 ], [ 0, %72 ], [ %108, %157 ]
  %95 = and i64 %73, 1
  %96 = icmp eq i64 %95, 0
  br i1 %96, label %97, label %104

97:                                               ; preds = %93
  %98 = add nuw nsw i64 %94, %84
  %99 = mul nsw i64 %98, %13
  %100 = add nuw nsw i64 %99, %18
  %101 = getelementptr inbounds float, float* %0, i64 %100
  store float 0.000000e+00, float* %101, align 4, !tbaa !8
  %102 = getelementptr inbounds float, float* %1, i64 %100
  store float 0.000000e+00, float* %102, align 4, !tbaa !8
  %103 = or i64 %94, 1
  br label %104

104:                                              ; preds = %97, %93
  %105 = phi i64 [ %103, %97 ], [ %94, %93 ]
  %106 = icmp eq i64 %73, %94
  br i1 %106, label %68, label %162

107:                                              ; preds = %87
  %108 = and i64 %80, -8
  %109 = and i64 %78, 1
  %110 = icmp eq i64 %77, 0
  br i1 %110, label %142, label %111

111:                                              ; preds = %107
  %112 = sub nuw nsw i64 %78, %109
  br label %113

113:                                              ; preds = %113, %111
  %114 = phi i64 [ 0, %111 ], [ %139, %113 ]
  %115 = phi i64 [ %112, %111 ], [ %140, %113 ]
  %116 = add nuw nsw i64 %114, %84
  %117 = mul nsw i64 %116, %13
  %118 = add nuw nsw i64 %117, %18
  %119 = getelementptr inbounds float, float* %0, i64 %118
  %120 = bitcast float* %119 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %120, align 4, !tbaa !8, !alias.scope !51, !noalias !54
  %121 = getelementptr inbounds float, float* %119, i64 4
  %122 = bitcast float* %121 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %122, align 4, !tbaa !8, !alias.scope !51, !noalias !54
  %123 = getelementptr inbounds float, float* %1, i64 %118
  %124 = bitcast float* %123 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %124, align 4, !tbaa !8, !alias.scope !54
  %125 = getelementptr inbounds float, float* %123, i64 4
  %126 = bitcast float* %125 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %126, align 4, !tbaa !8, !alias.scope !54
  %127 = or i64 %114, 8
  %128 = add nuw nsw i64 %127, %84
  %129 = mul nsw i64 %128, %13
  %130 = add nuw nsw i64 %129, %18
  %131 = getelementptr inbounds float, float* %0, i64 %130
  %132 = bitcast float* %131 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %132, align 4, !tbaa !8, !alias.scope !51, !noalias !54
  %133 = getelementptr inbounds float, float* %131, i64 4
  %134 = bitcast float* %133 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %134, align 4, !tbaa !8, !alias.scope !51, !noalias !54
  %135 = getelementptr inbounds float, float* %1, i64 %130
  %136 = bitcast float* %135 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %136, align 4, !tbaa !8, !alias.scope !54
  %137 = getelementptr inbounds float, float* %135, i64 4
  %138 = bitcast float* %137 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %138, align 4, !tbaa !8, !alias.scope !54
  %139 = add i64 %114, 16
  %140 = add i64 %115, -2
  %141 = icmp eq i64 %140, 0
  br i1 %141, label %142, label %113, !llvm.loop !56

142:                                              ; preds = %113, %107
  %143 = phi i64 [ 0, %107 ], [ %139, %113 ]
  %144 = icmp eq i64 %109, 0
  br i1 %144, label %157, label %145

145:                                              ; preds = %142
  %146 = add nuw nsw i64 %143, %84
  %147 = mul nsw i64 %146, %13
  %148 = add nuw nsw i64 %147, %18
  %149 = getelementptr inbounds float, float* %0, i64 %148
  %150 = bitcast float* %149 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %150, align 4, !tbaa !8, !alias.scope !51, !noalias !54
  %151 = getelementptr inbounds float, float* %149, i64 4
  %152 = bitcast float* %151 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %152, align 4, !tbaa !8, !alias.scope !51, !noalias !54
  %153 = getelementptr inbounds float, float* %1, i64 %148
  %154 = bitcast float* %153 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %154, align 4, !tbaa !8, !alias.scope !54
  %155 = getelementptr inbounds float, float* %153, i64 4
  %156 = bitcast float* %155 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %156, align 4, !tbaa !8, !alias.scope !54
  br label %157

157:                                              ; preds = %142, %145
  %158 = icmp eq i64 %80, %108
  br i1 %158, label %68, label %93

159:                                              ; preds = %68, %64, %17
  %160 = add nuw nsw i64 %18, 1
  %161 = icmp eq i64 %160, %13
  br i1 %161, label %67, label %17

162:                                              ; preds = %104, %162
  %163 = phi i64 [ %175, %162 ], [ %105, %104 ]
  %164 = add nuw nsw i64 %163, %84
  %165 = mul nsw i64 %164, %13
  %166 = add nuw nsw i64 %165, %18
  %167 = getelementptr inbounds float, float* %0, i64 %166
  store float 0.000000e+00, float* %167, align 4, !tbaa !8
  %168 = getelementptr inbounds float, float* %1, i64 %166
  store float 0.000000e+00, float* %168, align 4, !tbaa !8
  %169 = add nuw nsw i64 %163, 1
  %170 = add nuw nsw i64 %169, %84
  %171 = mul nsw i64 %170, %13
  %172 = add nuw nsw i64 %171, %18
  %173 = getelementptr inbounds float, float* %0, i64 %172
  store float 0.000000e+00, float* %173, align 4, !tbaa !8
  %174 = getelementptr inbounds float, float* %1, i64 %172
  store float 0.000000e+00, float* %174, align 4, !tbaa !8
  %175 = add nuw nsw i64 %163, 2
  %176 = icmp eq i64 %175, %74
  br i1 %176, label %68, label %162, !llvm.loop !57
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z38cpuWeightedRadialSphericalHarmonicsSumIdEvPT_S1_S1_S1_S1_S1_iii(double* %0, double* %1, double* %2, double* %3, double* %4, double* %5, i32 %6, i32 %7, i32 %8) local_unnamed_addr #0 {
  %10 = icmp sgt i32 %7, 0
  br i1 %10, label %11, label %73

11:                                               ; preds = %9
  %12 = icmp slt i32 %6, 0
  %13 = icmp sgt i32 %8, 0
  %14 = zext i32 %7 to i64
  %15 = sext i32 %8 to i64
  %16 = zext i32 %8 to i64
  %17 = icmp eq i32 %7, 1
  br label %18

18:                                               ; preds = %165, %11
  %19 = phi i64 [ 0, %11 ], [ %166, %165 ]
  br i1 %12, label %165, label %20

20:                                               ; preds = %18
  br i1 %13, label %21, label %78

21:                                               ; preds = %20
  %22 = trunc i64 %19 to i32
  br label %23

23:                                               ; preds = %70, %21
  %24 = phi i64 [ 1, %21 ], [ %71, %70 ]
  %25 = phi i32 [ 0, %21 ], [ %26, %70 ]
  %26 = add nuw i32 %25, 1
  %27 = mul nsw i32 %26, %25
  %28 = lshr i32 %27, 1
  %29 = mul i32 %25, %7
  %30 = add i32 %29, %22
  %31 = mul i32 %30, %8
  %32 = zext i32 %28 to i64
  br label %33

33:                                               ; preds = %41, %23
  %34 = phi i64 [ %42, %41 ], [ 0, %23 ]
  %35 = add nuw nsw i64 %34, %32
  %36 = mul nsw i64 %35, %14
  %37 = add nuw nsw i64 %36, %19
  %38 = getelementptr inbounds double, double* %0, i64 %37
  store double 0.000000e+00, double* %38, align 8, !tbaa !4
  %39 = getelementptr inbounds double, double* %1, i64 %37
  store double 0.000000e+00, double* %39, align 8, !tbaa !4
  %40 = mul nsw i64 %35, %15
  br label %44

41:                                               ; preds = %44
  %42 = add nuw nsw i64 %34, 1
  %43 = icmp eq i64 %42, %24
  br i1 %43, label %70, label %33

44:                                               ; preds = %44, %33
  %45 = phi i64 [ %68, %44 ], [ 0, %33 ]
  %46 = getelementptr inbounds double, double* %5, i64 %45
  %47 = load double, double* %46, align 8, !tbaa !4
  %48 = trunc i64 %45 to i32
  %49 = add i32 %31, %48
  %50 = sext i32 %49 to i64
  %51 = getelementptr inbounds double, double* %2, i64 %50
  %52 = load double, double* %51, align 8, !tbaa !4
  %53 = fmul fast double %52, %47
  %54 = add nsw i64 %45, %40
  %55 = getelementptr inbounds double, double* %3, i64 %54
  %56 = load double, double* %55, align 8, !tbaa !4
  %57 = fmul fast double %53, %56
  %58 = load double, double* %38, align 8, !tbaa !4
  %59 = fadd fast double %58, %57
  store double %59, double* %38, align 8, !tbaa !4
  %60 = load double, double* %46, align 8, !tbaa !4
  %61 = load double, double* %51, align 8, !tbaa !4
  %62 = fmul fast double %61, %60
  %63 = getelementptr inbounds double, double* %4, i64 %54
  %64 = load double, double* %63, align 8, !tbaa !4
  %65 = fmul fast double %62, %64
  %66 = load double, double* %39, align 8, !tbaa !4
  %67 = fadd fast double %66, %65
  store double %67, double* %39, align 8, !tbaa !4
  %68 = add nuw nsw i64 %45, 1
  %69 = icmp eq i64 %68, %16
  br i1 %69, label %41, label %44

70:                                               ; preds = %41
  %71 = add nuw nsw i64 %24, 1
  %72 = icmp eq i32 %25, %6
  br i1 %72, label %165, label %23

73:                                               ; preds = %165, %9
  ret void

74:                                               ; preds = %110, %168, %163
  %75 = add nuw nsw i64 %80, 1
  %76 = icmp eq i32 %81, %6
  %77 = add i64 %79, 1
  br i1 %76, label %165, label %78

78:                                               ; preds = %20, %74
  %79 = phi i64 [ %77, %74 ], [ 0, %20 ]
  %80 = phi i64 [ %75, %74 ], [ 1, %20 ]
  %81 = phi i32 [ %87, %74 ], [ 0, %20 ]
  %82 = add i64 %79, -3
  %83 = lshr i64 %82, 2
  %84 = add nuw nsw i64 %83, 1
  %85 = add i64 %79, 1
  %86 = add i64 %79, 1
  %87 = add nuw i32 %81, 1
  %88 = mul nsw i32 %87, %81
  %89 = lshr i32 %88, 1
  %90 = zext i32 %89 to i64
  %91 = icmp ugt i64 %86, 3
  %92 = and i1 %91, %17
  br i1 %92, label %93, label %99

93:                                               ; preds = %78
  %94 = getelementptr double, double* %1, i64 %85
  %95 = getelementptr double, double* %0, i64 %85
  %96 = icmp ugt double* %94, %0
  %97 = icmp ugt double* %95, %1
  %98 = and i1 %96, %97
  br i1 %98, label %99, label %113

99:                                               ; preds = %78, %163, %93
  %100 = phi i64 [ 0, %93 ], [ 0, %78 ], [ %114, %163 ]
  %101 = and i64 %79, 1
  %102 = icmp eq i64 %101, 0
  br i1 %102, label %103, label %110

103:                                              ; preds = %99
  %104 = add nuw nsw i64 %100, %90
  %105 = mul nsw i64 %104, %14
  %106 = add nuw nsw i64 %105, %19
  %107 = getelementptr inbounds double, double* %0, i64 %106
  store double 0.000000e+00, double* %107, align 8, !tbaa !4
  %108 = getelementptr inbounds double, double* %1, i64 %106
  store double 0.000000e+00, double* %108, align 8, !tbaa !4
  %109 = or i64 %100, 1
  br label %110

110:                                              ; preds = %103, %99
  %111 = phi i64 [ %109, %103 ], [ %100, %99 ]
  %112 = icmp eq i64 %79, %100
  br i1 %112, label %74, label %168

113:                                              ; preds = %93
  %114 = and i64 %86, -4
  %115 = and i64 %84, 1
  %116 = icmp eq i64 %83, 0
  br i1 %116, label %148, label %117

117:                                              ; preds = %113
  %118 = sub nuw nsw i64 %84, %115
  br label %119

119:                                              ; preds = %119, %117
  %120 = phi i64 [ 0, %117 ], [ %145, %119 ]
  %121 = phi i64 [ %118, %117 ], [ %146, %119 ]
  %122 = add nuw nsw i64 %120, %90
  %123 = mul nsw i64 %122, %14
  %124 = add nuw nsw i64 %123, %19
  %125 = getelementptr inbounds double, double* %0, i64 %124
  %126 = bitcast double* %125 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %126, align 8, !tbaa !4, !alias.scope !58, !noalias !61
  %127 = getelementptr inbounds double, double* %125, i64 2
  %128 = bitcast double* %127 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %128, align 8, !tbaa !4, !alias.scope !58, !noalias !61
  %129 = getelementptr inbounds double, double* %1, i64 %124
  %130 = bitcast double* %129 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %130, align 8, !tbaa !4, !alias.scope !61
  %131 = getelementptr inbounds double, double* %129, i64 2
  %132 = bitcast double* %131 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %132, align 8, !tbaa !4, !alias.scope !61
  %133 = or i64 %120, 4
  %134 = add nuw nsw i64 %133, %90
  %135 = mul nsw i64 %134, %14
  %136 = add nuw nsw i64 %135, %19
  %137 = getelementptr inbounds double, double* %0, i64 %136
  %138 = bitcast double* %137 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %138, align 8, !tbaa !4, !alias.scope !58, !noalias !61
  %139 = getelementptr inbounds double, double* %137, i64 2
  %140 = bitcast double* %139 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %140, align 8, !tbaa !4, !alias.scope !58, !noalias !61
  %141 = getelementptr inbounds double, double* %1, i64 %136
  %142 = bitcast double* %141 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %142, align 8, !tbaa !4, !alias.scope !61
  %143 = getelementptr inbounds double, double* %141, i64 2
  %144 = bitcast double* %143 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %144, align 8, !tbaa !4, !alias.scope !61
  %145 = add i64 %120, 8
  %146 = add i64 %121, -2
  %147 = icmp eq i64 %146, 0
  br i1 %147, label %148, label %119, !llvm.loop !63

148:                                              ; preds = %119, %113
  %149 = phi i64 [ 0, %113 ], [ %145, %119 ]
  %150 = icmp eq i64 %115, 0
  br i1 %150, label %163, label %151

151:                                              ; preds = %148
  %152 = add nuw nsw i64 %149, %90
  %153 = mul nsw i64 %152, %14
  %154 = add nuw nsw i64 %153, %19
  %155 = getelementptr inbounds double, double* %0, i64 %154
  %156 = bitcast double* %155 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %156, align 8, !tbaa !4, !alias.scope !58, !noalias !61
  %157 = getelementptr inbounds double, double* %155, i64 2
  %158 = bitcast double* %157 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %158, align 8, !tbaa !4, !alias.scope !58, !noalias !61
  %159 = getelementptr inbounds double, double* %1, i64 %154
  %160 = bitcast double* %159 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %160, align 8, !tbaa !4, !alias.scope !61
  %161 = getelementptr inbounds double, double* %159, i64 2
  %162 = bitcast double* %161 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %162, align 8, !tbaa !4, !alias.scope !61
  br label %163

163:                                              ; preds = %148, %151
  %164 = icmp eq i64 %86, %114
  br i1 %164, label %74, label %99

165:                                              ; preds = %74, %70, %18
  %166 = add nuw nsw i64 %19, 1
  %167 = icmp eq i64 %166, %14
  br i1 %167, label %73, label %18

168:                                              ; preds = %110, %168
  %169 = phi i64 [ %181, %168 ], [ %111, %110 ]
  %170 = add nuw nsw i64 %169, %90
  %171 = mul nsw i64 %170, %14
  %172 = add nuw nsw i64 %171, %19
  %173 = getelementptr inbounds double, double* %0, i64 %172
  store double 0.000000e+00, double* %173, align 8, !tbaa !4
  %174 = getelementptr inbounds double, double* %1, i64 %172
  store double 0.000000e+00, double* %174, align 8, !tbaa !4
  %175 = add nuw nsw i64 %169, 1
  %176 = add nuw nsw i64 %175, %90
  %177 = mul nsw i64 %176, %14
  %178 = add nuw nsw i64 %177, %19
  %179 = getelementptr inbounds double, double* %0, i64 %178
  store double 0.000000e+00, double* %179, align 8, !tbaa !4
  %180 = getelementptr inbounds double, double* %1, i64 %178
  store double 0.000000e+00, double* %180, align 8, !tbaa !4
  %181 = add nuw nsw i64 %169, 2
  %182 = icmp eq i64 %181, %80
  br i1 %182, label %74, label %168, !llvm.loop !64
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z38cpuWeightedRadialSphericalHarmonicsSumIfEvPT_S1_S1_S1_S1_S1_iii(float* %0, float* %1, float* %2, float* %3, float* %4, float* %5, i32 %6, i32 %7, i32 %8) local_unnamed_addr #0 {
  %10 = icmp sgt i32 %7, 0
  br i1 %10, label %11, label %73

11:                                               ; preds = %9
  %12 = icmp slt i32 %6, 0
  %13 = icmp sgt i32 %8, 0
  %14 = zext i32 %7 to i64
  %15 = sext i32 %8 to i64
  %16 = zext i32 %8 to i64
  %17 = icmp eq i32 %7, 1
  br label %18

18:                                               ; preds = %165, %11
  %19 = phi i64 [ 0, %11 ], [ %166, %165 ]
  br i1 %12, label %165, label %20

20:                                               ; preds = %18
  br i1 %13, label %21, label %78

21:                                               ; preds = %20
  %22 = trunc i64 %19 to i32
  br label %23

23:                                               ; preds = %70, %21
  %24 = phi i64 [ 1, %21 ], [ %71, %70 ]
  %25 = phi i32 [ 0, %21 ], [ %26, %70 ]
  %26 = add nuw i32 %25, 1
  %27 = mul nsw i32 %26, %25
  %28 = lshr i32 %27, 1
  %29 = mul i32 %25, %7
  %30 = add i32 %29, %22
  %31 = mul i32 %30, %8
  %32 = zext i32 %28 to i64
  br label %33

33:                                               ; preds = %41, %23
  %34 = phi i64 [ %42, %41 ], [ 0, %23 ]
  %35 = add nuw nsw i64 %34, %32
  %36 = mul nsw i64 %35, %14
  %37 = add nuw nsw i64 %36, %19
  %38 = getelementptr inbounds float, float* %0, i64 %37
  store float 0.000000e+00, float* %38, align 4, !tbaa !8
  %39 = getelementptr inbounds float, float* %1, i64 %37
  store float 0.000000e+00, float* %39, align 4, !tbaa !8
  %40 = mul nsw i64 %35, %15
  br label %44

41:                                               ; preds = %44
  %42 = add nuw nsw i64 %34, 1
  %43 = icmp eq i64 %42, %24
  br i1 %43, label %70, label %33

44:                                               ; preds = %44, %33
  %45 = phi i64 [ %68, %44 ], [ 0, %33 ]
  %46 = getelementptr inbounds float, float* %5, i64 %45
  %47 = load float, float* %46, align 4, !tbaa !8
  %48 = trunc i64 %45 to i32
  %49 = add i32 %31, %48
  %50 = sext i32 %49 to i64
  %51 = getelementptr inbounds float, float* %2, i64 %50
  %52 = load float, float* %51, align 4, !tbaa !8
  %53 = fmul fast float %52, %47
  %54 = add nsw i64 %45, %40
  %55 = getelementptr inbounds float, float* %3, i64 %54
  %56 = load float, float* %55, align 4, !tbaa !8
  %57 = fmul fast float %53, %56
  %58 = load float, float* %38, align 4, !tbaa !8
  %59 = fadd fast float %58, %57
  store float %59, float* %38, align 4, !tbaa !8
  %60 = load float, float* %46, align 4, !tbaa !8
  %61 = load float, float* %51, align 4, !tbaa !8
  %62 = fmul fast float %61, %60
  %63 = getelementptr inbounds float, float* %4, i64 %54
  %64 = load float, float* %63, align 4, !tbaa !8
  %65 = fmul fast float %62, %64
  %66 = load float, float* %39, align 4, !tbaa !8
  %67 = fadd fast float %66, %65
  store float %67, float* %39, align 4, !tbaa !8
  %68 = add nuw nsw i64 %45, 1
  %69 = icmp eq i64 %68, %16
  br i1 %69, label %41, label %44

70:                                               ; preds = %41
  %71 = add nuw nsw i64 %24, 1
  %72 = icmp eq i32 %25, %6
  br i1 %72, label %165, label %23

73:                                               ; preds = %165, %9
  ret void

74:                                               ; preds = %110, %168, %163
  %75 = add nuw nsw i64 %80, 1
  %76 = icmp eq i32 %81, %6
  %77 = add i64 %79, 1
  br i1 %76, label %165, label %78

78:                                               ; preds = %20, %74
  %79 = phi i64 [ %77, %74 ], [ 0, %20 ]
  %80 = phi i64 [ %75, %74 ], [ 1, %20 ]
  %81 = phi i32 [ %87, %74 ], [ 0, %20 ]
  %82 = add i64 %79, -7
  %83 = lshr i64 %82, 3
  %84 = add nuw nsw i64 %83, 1
  %85 = add i64 %79, 1
  %86 = add i64 %79, 1
  %87 = add nuw i32 %81, 1
  %88 = mul nsw i32 %87, %81
  %89 = lshr i32 %88, 1
  %90 = zext i32 %89 to i64
  %91 = icmp ugt i64 %86, 7
  %92 = and i1 %91, %17
  br i1 %92, label %93, label %99

93:                                               ; preds = %78
  %94 = getelementptr float, float* %1, i64 %85
  %95 = getelementptr float, float* %0, i64 %85
  %96 = icmp ugt float* %94, %0
  %97 = icmp ugt float* %95, %1
  %98 = and i1 %96, %97
  br i1 %98, label %99, label %113

99:                                               ; preds = %78, %163, %93
  %100 = phi i64 [ 0, %93 ], [ 0, %78 ], [ %114, %163 ]
  %101 = and i64 %79, 1
  %102 = icmp eq i64 %101, 0
  br i1 %102, label %103, label %110

103:                                              ; preds = %99
  %104 = add nuw nsw i64 %100, %90
  %105 = mul nsw i64 %104, %14
  %106 = add nuw nsw i64 %105, %19
  %107 = getelementptr inbounds float, float* %0, i64 %106
  store float 0.000000e+00, float* %107, align 4, !tbaa !8
  %108 = getelementptr inbounds float, float* %1, i64 %106
  store float 0.000000e+00, float* %108, align 4, !tbaa !8
  %109 = or i64 %100, 1
  br label %110

110:                                              ; preds = %103, %99
  %111 = phi i64 [ %109, %103 ], [ %100, %99 ]
  %112 = icmp eq i64 %79, %100
  br i1 %112, label %74, label %168

113:                                              ; preds = %93
  %114 = and i64 %86, -8
  %115 = and i64 %84, 1
  %116 = icmp eq i64 %83, 0
  br i1 %116, label %148, label %117

117:                                              ; preds = %113
  %118 = sub nuw nsw i64 %84, %115
  br label %119

119:                                              ; preds = %119, %117
  %120 = phi i64 [ 0, %117 ], [ %145, %119 ]
  %121 = phi i64 [ %118, %117 ], [ %146, %119 ]
  %122 = add nuw nsw i64 %120, %90
  %123 = mul nsw i64 %122, %14
  %124 = add nuw nsw i64 %123, %19
  %125 = getelementptr inbounds float, float* %0, i64 %124
  %126 = bitcast float* %125 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %126, align 4, !tbaa !8, !alias.scope !65, !noalias !68
  %127 = getelementptr inbounds float, float* %125, i64 4
  %128 = bitcast float* %127 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %128, align 4, !tbaa !8, !alias.scope !65, !noalias !68
  %129 = getelementptr inbounds float, float* %1, i64 %124
  %130 = bitcast float* %129 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %130, align 4, !tbaa !8, !alias.scope !68
  %131 = getelementptr inbounds float, float* %129, i64 4
  %132 = bitcast float* %131 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %132, align 4, !tbaa !8, !alias.scope !68
  %133 = or i64 %120, 8
  %134 = add nuw nsw i64 %133, %90
  %135 = mul nsw i64 %134, %14
  %136 = add nuw nsw i64 %135, %19
  %137 = getelementptr inbounds float, float* %0, i64 %136
  %138 = bitcast float* %137 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %138, align 4, !tbaa !8, !alias.scope !65, !noalias !68
  %139 = getelementptr inbounds float, float* %137, i64 4
  %140 = bitcast float* %139 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %140, align 4, !tbaa !8, !alias.scope !65, !noalias !68
  %141 = getelementptr inbounds float, float* %1, i64 %136
  %142 = bitcast float* %141 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %142, align 4, !tbaa !8, !alias.scope !68
  %143 = getelementptr inbounds float, float* %141, i64 4
  %144 = bitcast float* %143 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %144, align 4, !tbaa !8, !alias.scope !68
  %145 = add i64 %120, 16
  %146 = add i64 %121, -2
  %147 = icmp eq i64 %146, 0
  br i1 %147, label %148, label %119, !llvm.loop !70

148:                                              ; preds = %119, %113
  %149 = phi i64 [ 0, %113 ], [ %145, %119 ]
  %150 = icmp eq i64 %115, 0
  br i1 %150, label %163, label %151

151:                                              ; preds = %148
  %152 = add nuw nsw i64 %149, %90
  %153 = mul nsw i64 %152, %14
  %154 = add nuw nsw i64 %153, %19
  %155 = getelementptr inbounds float, float* %0, i64 %154
  %156 = bitcast float* %155 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %156, align 4, !tbaa !8, !alias.scope !65, !noalias !68
  %157 = getelementptr inbounds float, float* %155, i64 4
  %158 = bitcast float* %157 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %158, align 4, !tbaa !8, !alias.scope !65, !noalias !68
  %159 = getelementptr inbounds float, float* %1, i64 %154
  %160 = bitcast float* %159 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %160, align 4, !tbaa !8, !alias.scope !68
  %161 = getelementptr inbounds float, float* %159, i64 4
  %162 = bitcast float* %161 to <4 x float>*
  store <4 x float> zeroinitializer, <4 x float>* %162, align 4, !tbaa !8, !alias.scope !68
  br label %163

163:                                              ; preds = %148, %151
  %164 = icmp eq i64 %86, %114
  br i1 %164, label %74, label %99

165:                                              ; preds = %74, %70, %18
  %166 = add nuw nsw i64 %19, 1
  %167 = icmp eq i64 %166, %14
  br i1 %167, label %73, label %18

168:                                              ; preds = %110, %168
  %169 = phi i64 [ %181, %168 ], [ %111, %110 ]
  %170 = add nuw nsw i64 %169, %90
  %171 = mul nsw i64 %170, %14
  %172 = add nuw nsw i64 %171, %19
  %173 = getelementptr inbounds float, float* %0, i64 %172
  store float 0.000000e+00, float* %173, align 4, !tbaa !8
  %174 = getelementptr inbounds float, float* %1, i64 %172
  store float 0.000000e+00, float* %174, align 4, !tbaa !8
  %175 = add nuw nsw i64 %169, 1
  %176 = add nuw nsw i64 %175, %90
  %177 = mul nsw i64 %176, %14
  %178 = add nuw nsw i64 %177, %19
  %179 = getelementptr inbounds float, float* %0, i64 %178
  store float 0.000000e+00, float* %179, align 4, !tbaa !8
  %180 = getelementptr inbounds float, float* %1, i64 %178
  store float 0.000000e+00, float* %180, align 4, !tbaa !8
  %181 = add nuw nsw i64 %169, 2
  %182 = icmp eq i64 %181, %80
  br i1 %182, label %74, label %168, !llvm.loop !71
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z32cpuRadialSphericalHarmonicsDerivIdEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_iii(double* %0, double* %1, double* %2, double* %3, double* %4, double* %5, double* %6, double* %7, double* %8, double* %9, double* %10, double* %11, double* %12, double* %13, double* %14, double* %15, double* %16, i32 %17, i32 %18, i32 %19) local_unnamed_addr #0 {
  %21 = icmp sgt i32 %19, 0
  br i1 %21, label %22, label %30

22:                                               ; preds = %20
  %23 = icmp slt i32 %18, 1
  %24 = icmp slt i32 %17, 0
  %25 = mul i32 %19, %18
  %26 = zext i32 %19 to i64
  %27 = add i32 %17, 1
  %28 = zext i32 %27 to i64
  %29 = or i1 %23, %24
  br label %31

30:                                               ; preds = %65, %20
  ret void

31:                                               ; preds = %65, %22
  %32 = phi i64 [ 0, %22 ], [ %66, %65 ]
  %33 = getelementptr inbounds double, double* %6, i64 %32
  %34 = load double, double* %33, align 8, !tbaa !4
  %35 = fmul fast double %34, %34
  %36 = getelementptr inbounds double, double* %7, i64 %32
  %37 = load double, double* %36, align 8, !tbaa !4
  %38 = fmul fast double %37, %37
  %39 = fadd fast double %38, %35
  %40 = getelementptr inbounds double, double* %8, i64 %32
  %41 = load double, double* %40, align 8, !tbaa !4
  %42 = fmul fast double %41, %41
  %43 = fadd fast double %39, %42
  %44 = tail call fast double @llvm.sqrt.f64(double %43)
  %45 = tail call fast double @llvm.sqrt.f64(double %39)
  %46 = fmul fast double %45, %43
  %47 = fmul fast double %41, %34
  %48 = fdiv fast double %47, %46
  %49 = fmul fast double %41, %37
  %50 = fdiv fast double %49, %46
  %51 = fneg fast double %45
  %52 = fdiv fast double %51, %43
  %53 = fneg fast double %37
  %54 = fdiv fast double %53, %39
  %55 = fdiv fast double %34, %39
  br i1 %29, label %65, label %56

56:                                               ; preds = %31
  %57 = trunc i64 %32 to i32
  %58 = fdiv fast double 1.000000e+00, %44
  %59 = fdiv fast double 1.000000e+00, %44
  %60 = fdiv fast double 1.000000e+00, %44
  br label %61

61:                                               ; preds = %71, %56
  %62 = phi i32 [ %72, %71 ], [ 0, %56 ]
  %63 = mul nsw i32 %62, %19
  %64 = add i32 %63, %57
  br label %74

65:                                               ; preds = %71, %31
  %66 = add nuw nsw i64 %32, 1
  %67 = icmp eq i64 %66, %26
  br i1 %67, label %30, label %31

68:                                               ; preds = %96
  %69 = add nuw nsw i64 %76, 1
  %70 = icmp eq i64 %90, %28
  br i1 %70, label %71, label %74

71:                                               ; preds = %68
  %72 = add nuw nsw i32 %62, 1
  %73 = icmp eq i32 %72, %18
  br i1 %73, label %65, label %61

74:                                               ; preds = %68, %61
  %75 = phi i64 [ %90, %68 ], [ 0, %61 ]
  %76 = phi i64 [ %69, %68 ], [ 1, %61 ]
  %77 = phi i32 [ %91, %68 ], [ 0, %61 ]
  %78 = trunc i64 %75 to i32
  %79 = mul i32 %25, %78
  %80 = add i32 %64, %79
  %81 = sext i32 %80 to i64
  %82 = getelementptr inbounds double, double* %12, i64 %81
  %83 = load double, double* %82, align 8, !tbaa !4
  %84 = fmul fast double %83, %34
  %85 = fmul fast double %84, %58
  %86 = fmul fast double %83, %37
  %87 = fmul fast double %86, %59
  %88 = fmul fast double %83, %41
  %89 = fmul fast double %88, %60
  %90 = add nuw nsw i64 %75, 1
  %91 = add nuw nsw i32 %77, 1
  %92 = mul nsw i32 %91, %78
  %93 = lshr i32 %92, 1
  %94 = getelementptr inbounds double, double* %9, i64 %81
  %95 = zext i32 %93 to i64
  br label %96

96:                                               ; preds = %96, %74
  %97 = phi i64 [ %165, %96 ], [ 0, %74 ]
  %98 = add nuw nsw i64 %97, %95
  %99 = mul nsw i64 %98, %26
  %100 = trunc i64 %99 to i32
  %101 = mul i32 %100, %18
  %102 = add i32 %64, %101
  %103 = add nuw nsw i64 %99, %32
  %104 = getelementptr inbounds double, double* %13, i64 %103
  %105 = load double, double* %104, align 8, !tbaa !4
  %106 = fmul fast double %105, %48
  %107 = getelementptr inbounds double, double* %15, i64 %103
  %108 = load double, double* %107, align 8, !tbaa !4
  %109 = fmul fast double %108, %54
  %110 = fadd fast double %109, %106
  %111 = fmul fast double %105, %50
  %112 = fmul fast double %108, %55
  %113 = fadd fast double %112, %111
  %114 = fmul fast double %105, %52
  %115 = getelementptr inbounds double, double* %14, i64 %103
  %116 = load double, double* %115, align 8, !tbaa !4
  %117 = fmul fast double %116, %48
  %118 = getelementptr inbounds double, double* %16, i64 %103
  %119 = load double, double* %118, align 8, !tbaa !4
  %120 = fmul fast double %119, %54
  %121 = fadd fast double %120, %117
  %122 = fmul fast double %116, %50
  %123 = fmul fast double %119, %55
  %124 = fadd fast double %123, %122
  %125 = fmul fast double %116, %52
  %126 = getelementptr inbounds double, double* %10, i64 %103
  %127 = load double, double* %126, align 8, !tbaa !4
  %128 = fmul fast double %127, %85
  %129 = load double, double* %94, align 8, !tbaa !4
  %130 = fmul fast double %129, %110
  %131 = fadd fast double %130, %128
  %132 = sext i32 %102 to i64
  %133 = getelementptr inbounds double, double* %0, i64 %132
  store double %131, double* %133, align 8, !tbaa !4
  %134 = load double, double* %126, align 8, !tbaa !4
  %135 = fmul fast double %134, %87
  %136 = load double, double* %94, align 8, !tbaa !4
  %137 = fmul fast double %136, %113
  %138 = fadd fast double %137, %135
  %139 = getelementptr inbounds double, double* %2, i64 %132
  store double %138, double* %139, align 8, !tbaa !4
  %140 = load double, double* %126, align 8, !tbaa !4
  %141 = fmul fast double %140, %89
  %142 = load double, double* %94, align 8, !tbaa !4
  %143 = fmul fast double %114, %142
  %144 = fadd fast double %143, %141
  %145 = getelementptr inbounds double, double* %4, i64 %132
  store double %144, double* %145, align 8, !tbaa !4
  %146 = getelementptr inbounds double, double* %11, i64 %103
  %147 = load double, double* %146, align 8, !tbaa !4
  %148 = fmul fast double %147, %85
  %149 = load double, double* %94, align 8, !tbaa !4
  %150 = fmul fast double %149, %121
  %151 = fadd fast double %150, %148
  %152 = getelementptr inbounds double, double* %1, i64 %132
  store double %151, double* %152, align 8, !tbaa !4
  %153 = load double, double* %146, align 8, !tbaa !4
  %154 = fmul fast double %153, %87
  %155 = load double, double* %94, align 8, !tbaa !4
  %156 = fmul fast double %155, %124
  %157 = fadd fast double %156, %154
  %158 = getelementptr inbounds double, double* %3, i64 %132
  store double %157, double* %158, align 8, !tbaa !4
  %159 = load double, double* %146, align 8, !tbaa !4
  %160 = fmul fast double %159, %89
  %161 = load double, double* %94, align 8, !tbaa !4
  %162 = fmul fast double %125, %161
  %163 = fadd fast double %162, %160
  %164 = getelementptr inbounds double, double* %5, i64 %132
  store double %163, double* %164, align 8, !tbaa !4
  %165 = add nuw nsw i64 %97, 1
  %166 = icmp eq i64 %165, %76
  br i1 %166, label %68, label %96
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z32cpuRadialSphericalHarmonicsDerivIfEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_iii(float* %0, float* %1, float* %2, float* %3, float* %4, float* %5, float* %6, float* %7, float* %8, float* %9, float* %10, float* %11, float* %12, float* %13, float* %14, float* %15, float* %16, i32 %17, i32 %18, i32 %19) local_unnamed_addr #0 {
  %21 = icmp sgt i32 %19, 0
  br i1 %21, label %22, label %30

22:                                               ; preds = %20
  %23 = icmp slt i32 %18, 1
  %24 = icmp slt i32 %17, 0
  %25 = mul i32 %19, %18
  %26 = zext i32 %19 to i64
  %27 = add i32 %17, 1
  %28 = zext i32 %27 to i64
  %29 = or i1 %23, %24
  br label %31

30:                                               ; preds = %65, %20
  ret void

31:                                               ; preds = %65, %22
  %32 = phi i64 [ 0, %22 ], [ %66, %65 ]
  %33 = getelementptr inbounds float, float* %6, i64 %32
  %34 = load float, float* %33, align 4, !tbaa !8
  %35 = fmul fast float %34, %34
  %36 = getelementptr inbounds float, float* %7, i64 %32
  %37 = load float, float* %36, align 4, !tbaa !8
  %38 = fmul fast float %37, %37
  %39 = fadd fast float %38, %35
  %40 = getelementptr inbounds float, float* %8, i64 %32
  %41 = load float, float* %40, align 4, !tbaa !8
  %42 = fmul fast float %41, %41
  %43 = fadd fast float %39, %42
  %44 = tail call fast float @llvm.sqrt.f32(float %43) #8
  %45 = tail call fast float @llvm.sqrt.f32(float %39) #8
  %46 = fmul fast float %45, %43
  %47 = fmul fast float %41, %34
  %48 = fdiv fast float %47, %46
  %49 = fmul fast float %41, %37
  %50 = fdiv fast float %49, %46
  %51 = fneg fast float %45
  %52 = fdiv fast float %51, %43
  %53 = fneg fast float %37
  %54 = fdiv fast float %53, %39
  %55 = fdiv fast float %34, %39
  br i1 %29, label %65, label %56

56:                                               ; preds = %31
  %57 = trunc i64 %32 to i32
  %58 = fdiv fast float 1.000000e+00, %44
  %59 = fdiv fast float 1.000000e+00, %44
  %60 = fdiv fast float 1.000000e+00, %44
  br label %61

61:                                               ; preds = %71, %56
  %62 = phi i32 [ %72, %71 ], [ 0, %56 ]
  %63 = mul nsw i32 %62, %19
  %64 = add i32 %63, %57
  br label %74

65:                                               ; preds = %71, %31
  %66 = add nuw nsw i64 %32, 1
  %67 = icmp eq i64 %66, %26
  br i1 %67, label %30, label %31

68:                                               ; preds = %96
  %69 = add nuw nsw i64 %76, 1
  %70 = icmp eq i64 %90, %28
  br i1 %70, label %71, label %74

71:                                               ; preds = %68
  %72 = add nuw nsw i32 %62, 1
  %73 = icmp eq i32 %72, %18
  br i1 %73, label %65, label %61

74:                                               ; preds = %68, %61
  %75 = phi i64 [ %90, %68 ], [ 0, %61 ]
  %76 = phi i64 [ %69, %68 ], [ 1, %61 ]
  %77 = phi i32 [ %91, %68 ], [ 0, %61 ]
  %78 = trunc i64 %75 to i32
  %79 = mul i32 %25, %78
  %80 = add i32 %64, %79
  %81 = sext i32 %80 to i64
  %82 = getelementptr inbounds float, float* %12, i64 %81
  %83 = load float, float* %82, align 4, !tbaa !8
  %84 = fmul fast float %83, %34
  %85 = fmul fast float %84, %58
  %86 = fmul fast float %83, %37
  %87 = fmul fast float %86, %59
  %88 = fmul fast float %83, %41
  %89 = fmul fast float %88, %60
  %90 = add nuw nsw i64 %75, 1
  %91 = add nuw nsw i32 %77, 1
  %92 = mul nsw i32 %91, %78
  %93 = lshr i32 %92, 1
  %94 = getelementptr inbounds float, float* %9, i64 %81
  %95 = zext i32 %93 to i64
  br label %96

96:                                               ; preds = %96, %74
  %97 = phi i64 [ %165, %96 ], [ 0, %74 ]
  %98 = add nuw nsw i64 %97, %95
  %99 = mul nsw i64 %98, %26
  %100 = trunc i64 %99 to i32
  %101 = mul i32 %100, %18
  %102 = add i32 %64, %101
  %103 = add nuw nsw i64 %99, %32
  %104 = getelementptr inbounds float, float* %13, i64 %103
  %105 = load float, float* %104, align 4, !tbaa !8
  %106 = fmul fast float %105, %48
  %107 = getelementptr inbounds float, float* %15, i64 %103
  %108 = load float, float* %107, align 4, !tbaa !8
  %109 = fmul fast float %108, %54
  %110 = fadd fast float %109, %106
  %111 = fmul fast float %105, %50
  %112 = fmul fast float %108, %55
  %113 = fadd fast float %112, %111
  %114 = fmul fast float %105, %52
  %115 = getelementptr inbounds float, float* %14, i64 %103
  %116 = load float, float* %115, align 4, !tbaa !8
  %117 = fmul fast float %116, %48
  %118 = getelementptr inbounds float, float* %16, i64 %103
  %119 = load float, float* %118, align 4, !tbaa !8
  %120 = fmul fast float %119, %54
  %121 = fadd fast float %120, %117
  %122 = fmul fast float %116, %50
  %123 = fmul fast float %119, %55
  %124 = fadd fast float %123, %122
  %125 = fmul fast float %116, %52
  %126 = getelementptr inbounds float, float* %10, i64 %103
  %127 = load float, float* %126, align 4, !tbaa !8
  %128 = fmul fast float %127, %85
  %129 = load float, float* %94, align 4, !tbaa !8
  %130 = fmul fast float %129, %110
  %131 = fadd fast float %130, %128
  %132 = sext i32 %102 to i64
  %133 = getelementptr inbounds float, float* %0, i64 %132
  store float %131, float* %133, align 4, !tbaa !8
  %134 = load float, float* %126, align 4, !tbaa !8
  %135 = fmul fast float %134, %87
  %136 = load float, float* %94, align 4, !tbaa !8
  %137 = fmul fast float %136, %113
  %138 = fadd fast float %137, %135
  %139 = getelementptr inbounds float, float* %2, i64 %132
  store float %138, float* %139, align 4, !tbaa !8
  %140 = load float, float* %126, align 4, !tbaa !8
  %141 = fmul fast float %140, %89
  %142 = load float, float* %94, align 4, !tbaa !8
  %143 = fmul fast float %114, %142
  %144 = fadd fast float %143, %141
  %145 = getelementptr inbounds float, float* %4, i64 %132
  store float %144, float* %145, align 4, !tbaa !8
  %146 = getelementptr inbounds float, float* %11, i64 %103
  %147 = load float, float* %146, align 4, !tbaa !8
  %148 = fmul fast float %147, %85
  %149 = load float, float* %94, align 4, !tbaa !8
  %150 = fmul fast float %149, %121
  %151 = fadd fast float %150, %148
  %152 = getelementptr inbounds float, float* %1, i64 %132
  store float %151, float* %152, align 4, !tbaa !8
  %153 = load float, float* %146, align 4, !tbaa !8
  %154 = fmul fast float %153, %87
  %155 = load float, float* %94, align 4, !tbaa !8
  %156 = fmul fast float %155, %124
  %157 = fadd fast float %156, %154
  %158 = getelementptr inbounds float, float* %3, i64 %132
  store float %157, float* %158, align 4, !tbaa !8
  %159 = load float, float* %146, align 4, !tbaa !8
  %160 = fmul fast float %159, %89
  %161 = load float, float* %94, align 4, !tbaa !8
  %162 = fmul fast float %125, %161
  %163 = fadd fast float %162, %160
  %164 = getelementptr inbounds float, float* %5, i64 %132
  store float %163, float* %164, align 4, !tbaa !8
  %165 = add nuw nsw i64 %97, 1
  %166 = icmp eq i64 %165, %76
  br i1 %166, label %68, label %96
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z40cpuWeightedRadialSphericalHarmonicsDerivIdEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_iii(double* %0, double* %1, double* %2, double* %3, double* %4, double* %5, double* %6, double* %7, double* %8, double* %9, double* %10, double* %11, double* %12, double* %13, double* %14, double* %15, double* %16, double* %17, i32 %18, i32 %19, i32 %20) local_unnamed_addr #0 {
  %22 = icmp sgt i32 %20, 0
  br i1 %22, label %23, label %30

23:                                               ; preds = %21
  %24 = icmp sgt i32 %19, 0
  %25 = icmp slt i32 %18, 0
  %26 = mul i32 %20, %19
  %27 = zext i32 %20 to i64
  %28 = add i32 %18, 1
  %29 = zext i32 %28 to i64
  br label %31

30:                                               ; preds = %67, %21
  ret void

31:                                               ; preds = %67, %23
  %32 = phi i64 [ 0, %23 ], [ %68, %67 ]
  %33 = getelementptr inbounds double, double* %6, i64 %32
  %34 = load double, double* %33, align 8, !tbaa !4
  %35 = fmul fast double %34, %34
  %36 = getelementptr inbounds double, double* %7, i64 %32
  %37 = load double, double* %36, align 8, !tbaa !4
  %38 = fmul fast double %37, %37
  %39 = fadd fast double %38, %35
  %40 = getelementptr inbounds double, double* %8, i64 %32
  %41 = load double, double* %40, align 8, !tbaa !4
  %42 = fmul fast double %41, %41
  %43 = fadd fast double %39, %42
  %44 = tail call fast double @llvm.sqrt.f64(double %43)
  %45 = tail call fast double @llvm.sqrt.f64(double %39)
  %46 = fmul fast double %45, %43
  %47 = fmul fast double %41, %34
  %48 = fdiv fast double %47, %46
  %49 = fmul fast double %41, %37
  %50 = fdiv fast double %49, %46
  %51 = fneg fast double %45
  %52 = fneg fast double %37
  %53 = fdiv fast double %52, %39
  %54 = fdiv fast double %34, %39
  br i1 %24, label %55, label %67

55:                                               ; preds = %31
  %56 = getelementptr inbounds double, double* %17, i64 %32
  br i1 %25, label %67, label %57

57:                                               ; preds = %55
  %58 = trunc i64 %32 to i32
  %59 = fdiv fast double 1.000000e+00, %44
  %60 = fdiv fast double 1.000000e+00, %44
  %61 = fdiv fast double 1.000000e+00, %44
  %62 = fdiv fast double 1.000000e+00, %43
  br label %63

63:                                               ; preds = %73, %57
  %64 = phi i32 [ %74, %73 ], [ 0, %57 ]
  %65 = mul nsw i32 %64, %20
  %66 = add i32 %65, %58
  br label %76

67:                                               ; preds = %73, %55, %31
  %68 = add nuw nsw i64 %32, 1
  %69 = icmp eq i64 %68, %27
  br i1 %69, label %30, label %31

70:                                               ; preds = %104
  %71 = add nuw nsw i64 %78, 1
  %72 = icmp eq i64 %97, %29
  br i1 %72, label %73, label %76

73:                                               ; preds = %70
  %74 = add nuw nsw i32 %64, 1
  %75 = icmp eq i32 %74, %19
  br i1 %75, label %67, label %63

76:                                               ; preds = %70, %63
  %77 = phi i64 [ %97, %70 ], [ 0, %63 ]
  %78 = phi i64 [ %71, %70 ], [ 1, %63 ]
  %79 = phi i32 [ %98, %70 ], [ 0, %63 ]
  %80 = trunc i64 %77 to i32
  %81 = mul i32 %26, %80
  %82 = add i32 %66, %81
  %83 = load double, double* %56, align 8, !tbaa !4
  %84 = sext i32 %82 to i64
  %85 = getelementptr inbounds double, double* %12, i64 %84
  %86 = load double, double* %85, align 8, !tbaa !4
  %87 = fmul fast double %86, %83
  %88 = fmul fast double %87, %34
  %89 = fmul fast double %88, %59
  %90 = fmul fast double %87, %37
  %91 = fmul fast double %90, %60
  %92 = fmul fast double %87, %41
  %93 = fmul fast double %92, %61
  %94 = getelementptr inbounds double, double* %9, i64 %84
  %95 = load double, double* %94, align 8, !tbaa !4
  %96 = fmul fast double %95, %83
  %97 = add nuw nsw i64 %77, 1
  %98 = add nuw nsw i32 %79, 1
  %99 = mul nsw i32 %98, %80
  %100 = lshr i32 %99, 1
  %101 = fmul fast double %96, %51
  %102 = fmul fast double %101, %62
  %103 = zext i32 %100 to i64
  br label %104

104:                                              ; preds = %104, %76
  %105 = phi i64 [ %165, %104 ], [ 0, %76 ]
  %106 = add nuw nsw i64 %105, %103
  %107 = mul nsw i64 %106, %27
  %108 = trunc i64 %107 to i32
  %109 = mul i32 %108, %19
  %110 = add i32 %66, %109
  %111 = add nuw nsw i64 %107, %32
  %112 = getelementptr inbounds double, double* %13, i64 %111
  %113 = load double, double* %112, align 8, !tbaa !4
  %114 = fmul fast double %113, %48
  %115 = getelementptr inbounds double, double* %15, i64 %111
  %116 = load double, double* %115, align 8, !tbaa !4
  %117 = fmul fast double %116, %53
  %118 = fadd fast double %117, %114
  %119 = fmul fast double %113, %50
  %120 = fmul fast double %116, %54
  %121 = fadd fast double %120, %119
  %122 = getelementptr inbounds double, double* %14, i64 %111
  %123 = load double, double* %122, align 8, !tbaa !4
  %124 = fmul fast double %123, %48
  %125 = getelementptr inbounds double, double* %16, i64 %111
  %126 = load double, double* %125, align 8, !tbaa !4
  %127 = fmul fast double %126, %53
  %128 = fadd fast double %127, %124
  %129 = fmul fast double %123, %50
  %130 = fmul fast double %126, %54
  %131 = fadd fast double %130, %129
  %132 = getelementptr inbounds double, double* %10, i64 %111
  %133 = load double, double* %132, align 8, !tbaa !4
  %134 = fmul fast double %133, %89
  %135 = fmul fast double %118, %96
  %136 = fadd fast double %134, %135
  %137 = sext i32 %110 to i64
  %138 = getelementptr inbounds double, double* %0, i64 %137
  store double %136, double* %138, align 8, !tbaa !4
  %139 = load double, double* %132, align 8, !tbaa !4
  %140 = fmul fast double %139, %91
  %141 = fmul fast double %121, %96
  %142 = fadd fast double %140, %141
  %143 = getelementptr inbounds double, double* %2, i64 %137
  store double %142, double* %143, align 8, !tbaa !4
  %144 = load double, double* %132, align 8, !tbaa !4
  %145 = fmul fast double %144, %93
  %146 = fmul fast double %102, %113
  %147 = fadd fast double %145, %146
  %148 = getelementptr inbounds double, double* %4, i64 %137
  store double %147, double* %148, align 8, !tbaa !4
  %149 = getelementptr inbounds double, double* %11, i64 %111
  %150 = load double, double* %149, align 8, !tbaa !4
  %151 = fmul fast double %150, %89
  %152 = fmul fast double %128, %96
  %153 = fadd fast double %151, %152
  %154 = getelementptr inbounds double, double* %1, i64 %137
  store double %153, double* %154, align 8, !tbaa !4
  %155 = load double, double* %149, align 8, !tbaa !4
  %156 = fmul fast double %155, %91
  %157 = fmul fast double %131, %96
  %158 = fadd fast double %156, %157
  %159 = getelementptr inbounds double, double* %3, i64 %137
  store double %158, double* %159, align 8, !tbaa !4
  %160 = load double, double* %149, align 8, !tbaa !4
  %161 = fmul fast double %160, %93
  %162 = fmul fast double %102, %123
  %163 = fadd fast double %161, %162
  %164 = getelementptr inbounds double, double* %5, i64 %137
  store double %163, double* %164, align 8, !tbaa !4
  %165 = add nuw nsw i64 %105, 1
  %166 = icmp eq i64 %165, %78
  br i1 %166, label %70, label %104
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z40cpuWeightedRadialSphericalHarmonicsDerivIfEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_iii(float* %0, float* %1, float* %2, float* %3, float* %4, float* %5, float* %6, float* %7, float* %8, float* %9, float* %10, float* %11, float* %12, float* %13, float* %14, float* %15, float* %16, float* %17, i32 %18, i32 %19, i32 %20) local_unnamed_addr #0 {
  %22 = icmp sgt i32 %20, 0
  br i1 %22, label %23, label %30

23:                                               ; preds = %21
  %24 = icmp sgt i32 %19, 0
  %25 = icmp slt i32 %18, 0
  %26 = mul i32 %20, %19
  %27 = zext i32 %20 to i64
  %28 = add i32 %18, 1
  %29 = zext i32 %28 to i64
  br label %31

30:                                               ; preds = %67, %21
  ret void

31:                                               ; preds = %67, %23
  %32 = phi i64 [ 0, %23 ], [ %68, %67 ]
  %33 = getelementptr inbounds float, float* %6, i64 %32
  %34 = load float, float* %33, align 4, !tbaa !8
  %35 = fmul fast float %34, %34
  %36 = getelementptr inbounds float, float* %7, i64 %32
  %37 = load float, float* %36, align 4, !tbaa !8
  %38 = fmul fast float %37, %37
  %39 = fadd fast float %38, %35
  %40 = getelementptr inbounds float, float* %8, i64 %32
  %41 = load float, float* %40, align 4, !tbaa !8
  %42 = fmul fast float %41, %41
  %43 = fadd fast float %39, %42
  %44 = tail call fast float @llvm.sqrt.f32(float %43) #8
  %45 = tail call fast float @llvm.sqrt.f32(float %39) #8
  %46 = fmul fast float %45, %43
  %47 = fmul fast float %41, %34
  %48 = fdiv fast float %47, %46
  %49 = fmul fast float %41, %37
  %50 = fdiv fast float %49, %46
  %51 = fneg fast float %45
  %52 = fneg fast float %37
  %53 = fdiv fast float %52, %39
  %54 = fdiv fast float %34, %39
  br i1 %24, label %55, label %67

55:                                               ; preds = %31
  %56 = getelementptr inbounds float, float* %17, i64 %32
  br i1 %25, label %67, label %57

57:                                               ; preds = %55
  %58 = trunc i64 %32 to i32
  %59 = fdiv fast float 1.000000e+00, %44
  %60 = fdiv fast float 1.000000e+00, %44
  %61 = fdiv fast float 1.000000e+00, %44
  %62 = fdiv fast float 1.000000e+00, %43
  br label %63

63:                                               ; preds = %73, %57
  %64 = phi i32 [ %74, %73 ], [ 0, %57 ]
  %65 = mul nsw i32 %64, %20
  %66 = add i32 %65, %58
  br label %76

67:                                               ; preds = %73, %55, %31
  %68 = add nuw nsw i64 %32, 1
  %69 = icmp eq i64 %68, %27
  br i1 %69, label %30, label %31

70:                                               ; preds = %104
  %71 = add nuw nsw i64 %78, 1
  %72 = icmp eq i64 %97, %29
  br i1 %72, label %73, label %76

73:                                               ; preds = %70
  %74 = add nuw nsw i32 %64, 1
  %75 = icmp eq i32 %74, %19
  br i1 %75, label %67, label %63

76:                                               ; preds = %70, %63
  %77 = phi i64 [ %97, %70 ], [ 0, %63 ]
  %78 = phi i64 [ %71, %70 ], [ 1, %63 ]
  %79 = phi i32 [ %98, %70 ], [ 0, %63 ]
  %80 = trunc i64 %77 to i32
  %81 = mul i32 %26, %80
  %82 = add i32 %66, %81
  %83 = load float, float* %56, align 4, !tbaa !8
  %84 = sext i32 %82 to i64
  %85 = getelementptr inbounds float, float* %12, i64 %84
  %86 = load float, float* %85, align 4, !tbaa !8
  %87 = fmul fast float %86, %83
  %88 = fmul fast float %87, %34
  %89 = fmul fast float %88, %59
  %90 = fmul fast float %87, %37
  %91 = fmul fast float %90, %60
  %92 = fmul fast float %87, %41
  %93 = fmul fast float %92, %61
  %94 = getelementptr inbounds float, float* %9, i64 %84
  %95 = load float, float* %94, align 4, !tbaa !8
  %96 = fmul fast float %95, %83
  %97 = add nuw nsw i64 %77, 1
  %98 = add nuw nsw i32 %79, 1
  %99 = mul nsw i32 %98, %80
  %100 = lshr i32 %99, 1
  %101 = fmul fast float %96, %51
  %102 = fmul fast float %101, %62
  %103 = zext i32 %100 to i64
  br label %104

104:                                              ; preds = %104, %76
  %105 = phi i64 [ %165, %104 ], [ 0, %76 ]
  %106 = add nuw nsw i64 %105, %103
  %107 = mul nsw i64 %106, %27
  %108 = trunc i64 %107 to i32
  %109 = mul i32 %108, %19
  %110 = add i32 %66, %109
  %111 = add nuw nsw i64 %107, %32
  %112 = getelementptr inbounds float, float* %13, i64 %111
  %113 = load float, float* %112, align 4, !tbaa !8
  %114 = fmul fast float %113, %48
  %115 = getelementptr inbounds float, float* %15, i64 %111
  %116 = load float, float* %115, align 4, !tbaa !8
  %117 = fmul fast float %116, %53
  %118 = fadd fast float %117, %114
  %119 = fmul fast float %113, %50
  %120 = fmul fast float %116, %54
  %121 = fadd fast float %120, %119
  %122 = getelementptr inbounds float, float* %14, i64 %111
  %123 = load float, float* %122, align 4, !tbaa !8
  %124 = fmul fast float %123, %48
  %125 = getelementptr inbounds float, float* %16, i64 %111
  %126 = load float, float* %125, align 4, !tbaa !8
  %127 = fmul fast float %126, %53
  %128 = fadd fast float %127, %124
  %129 = fmul fast float %123, %50
  %130 = fmul fast float %126, %54
  %131 = fadd fast float %130, %129
  %132 = getelementptr inbounds float, float* %10, i64 %111
  %133 = load float, float* %132, align 4, !tbaa !8
  %134 = fmul fast float %133, %89
  %135 = fmul fast float %118, %96
  %136 = fadd fast float %134, %135
  %137 = sext i32 %110 to i64
  %138 = getelementptr inbounds float, float* %0, i64 %137
  store float %136, float* %138, align 4, !tbaa !8
  %139 = load float, float* %132, align 4, !tbaa !8
  %140 = fmul fast float %139, %91
  %141 = fmul fast float %121, %96
  %142 = fadd fast float %140, %141
  %143 = getelementptr inbounds float, float* %2, i64 %137
  store float %142, float* %143, align 4, !tbaa !8
  %144 = load float, float* %132, align 4, !tbaa !8
  %145 = fmul fast float %144, %93
  %146 = fmul fast float %102, %113
  %147 = fadd fast float %145, %146
  %148 = getelementptr inbounds float, float* %4, i64 %137
  store float %147, float* %148, align 4, !tbaa !8
  %149 = getelementptr inbounds float, float* %11, i64 %111
  %150 = load float, float* %149, align 4, !tbaa !8
  %151 = fmul fast float %150, %89
  %152 = fmul fast float %128, %96
  %153 = fadd fast float %151, %152
  %154 = getelementptr inbounds float, float* %1, i64 %137
  store float %153, float* %154, align 4, !tbaa !8
  %155 = load float, float* %149, align 4, !tbaa !8
  %156 = fmul fast float %155, %91
  %157 = fmul fast float %131, %96
  %158 = fadd fast float %156, %157
  %159 = getelementptr inbounds float, float* %3, i64 %137
  store float %158, float* %159, align 4, !tbaa !8
  %160 = load float, float* %149, align 4, !tbaa !8
  %161 = fmul fast float %160, %93
  %162 = fmul fast float %102, %123
  %163 = fadd fast float %161, %162
  %164 = getelementptr inbounds float, float* %5, i64 %137
  store float %163, float* %164, align 4, !tbaa !8
  %165 = add nuw nsw i64 %105, 1
  %166 = icmp eq i64 %165, %78
  br i1 %166, label %70, label %104
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z32cpuRadialSphericalHarmonicsPowerIdEvPT_S1_S1_Piii(double* %0, double* %1, double* %2, i32* %3, i32 %4, i32 %5) local_unnamed_addr #0 {
  %7 = add nsw i32 %5, 1
  %8 = mul nsw i32 %7, %5
  %9 = sdiv i32 %8, 2
  %10 = icmp sgt i32 %8, 1
  br i1 %10, label %11, label %20

11:                                               ; preds = %6
  %12 = icmp slt i32 %4, 0
  %13 = add i32 %4, 1
  br i1 %12, label %20, label %14

14:                                               ; preds = %11
  %15 = sext i32 %5 to i64
  %16 = sext i32 %13 to i64
  %17 = sext i32 %9 to i64
  %18 = zext i32 %9 to i64
  %19 = zext i32 %13 to i64
  br label %21

20:                                               ; preds = %34, %11, %6
  ret void

21:                                               ; preds = %34, %14
  %22 = phi i64 [ 0, %14 ], [ %35, %34 ]
  %23 = getelementptr inbounds i32, i32* %3, i64 %22
  %24 = load i32, i32* %23, align 4, !tbaa !72
  %25 = add nsw i64 %22, %17
  %26 = getelementptr inbounds i32, i32* %3, i64 %25
  %27 = load i32, i32* %26, align 4, !tbaa !72
  %28 = mul nsw i64 %22, %16
  %29 = sext i32 %27 to i64
  %30 = sext i32 %24 to i64
  br label %37

31:                                               ; preds = %67, %37
  %32 = add nuw nsw i64 %39, 1
  %33 = icmp eq i64 %41, %19
  br i1 %33, label %34, label %37

34:                                               ; preds = %31
  %35 = add nuw nsw i64 %22, 1
  %36 = icmp eq i64 %35, %18
  br i1 %36, label %20, label %21

37:                                               ; preds = %31, %21
  %38 = phi i64 [ %41, %31 ], [ 0, %21 ]
  %39 = phi i64 [ %32, %31 ], [ 1, %21 ]
  %40 = phi i32 [ %42, %31 ], [ 0, %21 ]
  %41 = add nuw nsw i64 %38, 1
  %42 = add nuw nsw i32 %40, 1
  %43 = trunc i64 %38 to i32
  %44 = mul nsw i32 %42, %43
  %45 = lshr i32 %44, 1
  %46 = mul nsw i32 %45, %5
  %47 = add nsw i32 %46, %27
  %48 = add nsw i32 %46, %24
  %49 = add nsw i64 %38, %28
  %50 = sext i32 %47 to i64
  %51 = getelementptr inbounds double, double* %1, i64 %50
  %52 = load double, double* %51, align 8, !tbaa !4
  %53 = sext i32 %48 to i64
  %54 = getelementptr inbounds double, double* %1, i64 %53
  %55 = load double, double* %54, align 8, !tbaa !4
  %56 = fmul fast double %55, %52
  %57 = getelementptr inbounds double, double* %2, i64 %50
  %58 = load double, double* %57, align 8, !tbaa !4
  %59 = getelementptr inbounds double, double* %2, i64 %53
  %60 = load double, double* %59, align 8, !tbaa !4
  %61 = fmul fast double %60, %58
  %62 = fadd fast double %61, %56
  %63 = getelementptr inbounds double, double* %0, i64 %49
  store double %62, double* %63, align 8, !tbaa !4
  %64 = icmp eq i64 %38, 0
  br i1 %64, label %31, label %65

65:                                               ; preds = %37
  %66 = zext i32 %45 to i64
  br label %67

67:                                               ; preds = %67, %65
  %68 = phi double [ %62, %65 ], [ %86, %67 ]
  %69 = phi i64 [ 1, %65 ], [ %87, %67 ]
  %70 = add nuw nsw i64 %69, %66
  %71 = mul nsw i64 %70, %15
  %72 = add nsw i64 %71, %29
  %73 = add nsw i64 %71, %30
  %74 = getelementptr inbounds double, double* %1, i64 %72
  %75 = load double, double* %74, align 8, !tbaa !4
  %76 = getelementptr inbounds double, double* %1, i64 %73
  %77 = load double, double* %76, align 8, !tbaa !4
  %78 = fmul fast double %77, %75
  %79 = getelementptr inbounds double, double* %2, i64 %72
  %80 = load double, double* %79, align 8, !tbaa !4
  %81 = getelementptr inbounds double, double* %2, i64 %73
  %82 = load double, double* %81, align 8, !tbaa !4
  %83 = fmul fast double %82, %80
  %84 = fadd fast double %83, %78
  %85 = fmul fast double %84, 2.000000e+00
  %86 = fadd fast double %85, %68
  store double %86, double* %63, align 8, !tbaa !4
  %87 = add nuw nsw i64 %69, 1
  %88 = icmp eq i64 %87, %39
  br i1 %88, label %31, label %67
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z32cpuRadialSphericalHarmonicsPowerIfEvPT_S1_S1_Piii(float* %0, float* %1, float* %2, i32* %3, i32 %4, i32 %5) local_unnamed_addr #0 {
  %7 = add nsw i32 %5, 1
  %8 = mul nsw i32 %7, %5
  %9 = sdiv i32 %8, 2
  %10 = icmp sgt i32 %8, 1
  br i1 %10, label %11, label %20

11:                                               ; preds = %6
  %12 = icmp slt i32 %4, 0
  %13 = add i32 %4, 1
  br i1 %12, label %20, label %14

14:                                               ; preds = %11
  %15 = sext i32 %5 to i64
  %16 = sext i32 %13 to i64
  %17 = sext i32 %9 to i64
  %18 = zext i32 %9 to i64
  %19 = zext i32 %13 to i64
  br label %21

20:                                               ; preds = %34, %11, %6
  ret void

21:                                               ; preds = %34, %14
  %22 = phi i64 [ 0, %14 ], [ %35, %34 ]
  %23 = getelementptr inbounds i32, i32* %3, i64 %22
  %24 = load i32, i32* %23, align 4, !tbaa !72
  %25 = add nsw i64 %22, %17
  %26 = getelementptr inbounds i32, i32* %3, i64 %25
  %27 = load i32, i32* %26, align 4, !tbaa !72
  %28 = mul nsw i64 %22, %16
  %29 = sext i32 %27 to i64
  %30 = sext i32 %24 to i64
  br label %37

31:                                               ; preds = %67, %37
  %32 = add nuw nsw i64 %39, 1
  %33 = icmp eq i64 %41, %19
  br i1 %33, label %34, label %37

34:                                               ; preds = %31
  %35 = add nuw nsw i64 %22, 1
  %36 = icmp eq i64 %35, %18
  br i1 %36, label %20, label %21

37:                                               ; preds = %31, %21
  %38 = phi i64 [ %41, %31 ], [ 0, %21 ]
  %39 = phi i64 [ %32, %31 ], [ 1, %21 ]
  %40 = phi i32 [ %42, %31 ], [ 0, %21 ]
  %41 = add nuw nsw i64 %38, 1
  %42 = add nuw nsw i32 %40, 1
  %43 = trunc i64 %38 to i32
  %44 = mul nsw i32 %42, %43
  %45 = lshr i32 %44, 1
  %46 = mul nsw i32 %45, %5
  %47 = add nsw i32 %46, %27
  %48 = add nsw i32 %46, %24
  %49 = add nsw i64 %38, %28
  %50 = sext i32 %47 to i64
  %51 = getelementptr inbounds float, float* %1, i64 %50
  %52 = load float, float* %51, align 4, !tbaa !8
  %53 = sext i32 %48 to i64
  %54 = getelementptr inbounds float, float* %1, i64 %53
  %55 = load float, float* %54, align 4, !tbaa !8
  %56 = fmul fast float %55, %52
  %57 = getelementptr inbounds float, float* %2, i64 %50
  %58 = load float, float* %57, align 4, !tbaa !8
  %59 = getelementptr inbounds float, float* %2, i64 %53
  %60 = load float, float* %59, align 4, !tbaa !8
  %61 = fmul fast float %60, %58
  %62 = fadd fast float %61, %56
  %63 = getelementptr inbounds float, float* %0, i64 %49
  store float %62, float* %63, align 4, !tbaa !8
  %64 = icmp eq i64 %38, 0
  br i1 %64, label %31, label %65

65:                                               ; preds = %37
  %66 = zext i32 %45 to i64
  br label %67

67:                                               ; preds = %67, %65
  %68 = phi float [ %62, %65 ], [ %89, %67 ]
  %69 = phi i64 [ 1, %65 ], [ %90, %67 ]
  %70 = add nuw nsw i64 %69, %66
  %71 = mul nsw i64 %70, %15
  %72 = add nsw i64 %71, %29
  %73 = add nsw i64 %71, %30
  %74 = getelementptr inbounds float, float* %1, i64 %72
  %75 = load float, float* %74, align 4, !tbaa !8
  %76 = getelementptr inbounds float, float* %1, i64 %73
  %77 = load float, float* %76, align 4, !tbaa !8
  %78 = fmul fast float %77, %75
  %79 = getelementptr inbounds float, float* %2, i64 %72
  %80 = load float, float* %79, align 4, !tbaa !8
  %81 = getelementptr inbounds float, float* %2, i64 %73
  %82 = load float, float* %81, align 4, !tbaa !8
  %83 = fmul fast float %82, %80
  %84 = fadd fast float %83, %78
  %85 = fpext float %84 to double
  %86 = fmul fast double %85, 2.000000e+00
  %87 = fpext float %68 to double
  %88 = fadd fast double %86, %87
  %89 = fptrunc double %88 to float
  store float %89, float* %63, align 4, !tbaa !8
  %90 = add nuw nsw i64 %69, 1
  %91 = icmp eq i64 %90, %39
  br i1 %91, label %31, label %67
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z37cpuRadialSphericalHarmonicsPowerDerivIdEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_Piiii(double* %0, double* %1, double* %2, double* %3, double* %4, double* %5, double* %6, double* %7, double* %8, double* %9, double* %10, i32* %11, i32 %12, i32 %13, i32 %14) local_unnamed_addr #0 {
  %16 = add nsw i32 %13, 1
  %17 = mul nsw i32 %16, %13
  %18 = sdiv i32 %17, 2
  %19 = icmp sgt i32 %14, 0
  br i1 %19, label %20, label %35

20:                                               ; preds = %15
  %21 = icmp sgt i32 %17, 1
  %22 = icmp slt i32 %12, 0
  %23 = add i32 %12, 1
  %24 = mul i32 %14, %13
  %25 = sext i32 %13 to i64
  %26 = zext i32 %14 to i64
  %27 = sext i32 %23 to i64
  %28 = sext i32 %18 to i64
  %29 = zext i32 %18 to i64
  %30 = zext i32 %23 to i64
  br label %31

31:                                               ; preds = %36, %20
  %32 = phi i64 [ 0, %20 ], [ %37, %36 ]
  br i1 %21, label %33, label %36

33:                                               ; preds = %31
  %34 = trunc i64 %32 to i32
  br label %39

35:                                               ; preds = %36, %15
  ret void

36:                                               ; preds = %57, %31
  %37 = add nuw nsw i64 %32, 1
  %38 = icmp eq i64 %37, %26
  br i1 %38, label %35, label %31

39:                                               ; preds = %57, %33
  %40 = phi i64 [ 0, %33 ], [ %58, %57 ]
  %41 = getelementptr inbounds i32, i32* %11, i64 %40
  %42 = load i32, i32* %41, align 4, !tbaa !72
  %43 = add nsw i64 %40, %28
  %44 = getelementptr inbounds i32, i32* %11, i64 %43
  %45 = load i32, i32* %44, align 4, !tbaa !72
  br i1 %22, label %57, label %46

46:                                               ; preds = %39
  %47 = mul nsw i64 %40, %27
  %48 = mul nsw i32 %45, %14
  %49 = add i32 %48, %34
  %50 = mul nsw i32 %42, %14
  %51 = add i32 %50, %34
  %52 = sext i32 %45 to i64
  %53 = sext i32 %42 to i64
  br label %60

54:                                               ; preds = %149, %60
  %55 = add nuw nsw i64 %62, 1
  %56 = icmp eq i64 %64, %30
  br i1 %56, label %57, label %60

57:                                               ; preds = %54, %39
  %58 = add nuw nsw i64 %40, 1
  %59 = icmp eq i64 %58, %29
  br i1 %59, label %36, label %39

60:                                               ; preds = %54, %46
  %61 = phi i64 [ %64, %54 ], [ 0, %46 ]
  %62 = phi i64 [ %55, %54 ], [ 1, %46 ]
  %63 = phi i32 [ %65, %54 ], [ 0, %46 ]
  %64 = add nuw nsw i64 %61, 1
  %65 = add nuw nsw i32 %63, 1
  %66 = trunc i64 %61 to i32
  %67 = mul nsw i32 %65, %66
  %68 = lshr i32 %67, 1
  %69 = mul nsw i32 %68, %13
  %70 = add nsw i32 %69, %45
  %71 = add nsw i32 %69, %42
  %72 = add nsw i64 %61, %47
  %73 = mul i32 %24, %68
  %74 = add i32 %49, %73
  %75 = add i32 %51, %73
  %76 = sext i32 %70 to i64
  %77 = getelementptr inbounds double, double* %3, i64 %76
  %78 = load double, double* %77, align 8, !tbaa !4
  %79 = sext i32 %75 to i64
  %80 = getelementptr inbounds double, double* %5, i64 %79
  %81 = load double, double* %80, align 8, !tbaa !4
  %82 = fmul fast double %81, %78
  %83 = sext i32 %74 to i64
  %84 = getelementptr inbounds double, double* %5, i64 %83
  %85 = load double, double* %84, align 8, !tbaa !4
  %86 = sext i32 %71 to i64
  %87 = getelementptr inbounds double, double* %3, i64 %86
  %88 = load double, double* %87, align 8, !tbaa !4
  %89 = fmul fast double %88, %85
  %90 = fadd fast double %89, %82
  %91 = getelementptr inbounds double, double* %4, i64 %76
  %92 = load double, double* %91, align 8, !tbaa !4
  %93 = getelementptr inbounds double, double* %6, i64 %79
  %94 = load double, double* %93, align 8, !tbaa !4
  %95 = fmul fast double %94, %92
  %96 = fadd fast double %90, %95
  %97 = getelementptr inbounds double, double* %6, i64 %83
  %98 = load double, double* %97, align 8, !tbaa !4
  %99 = getelementptr inbounds double, double* %4, i64 %86
  %100 = load double, double* %99, align 8, !tbaa !4
  %101 = fmul fast double %100, %98
  %102 = fadd fast double %96, %101
  %103 = mul nsw i64 %72, %26
  %104 = add nsw i64 %103, %32
  %105 = getelementptr inbounds double, double* %0, i64 %104
  store double %102, double* %105, align 8, !tbaa !4
  %106 = load double, double* %77, align 8, !tbaa !4
  %107 = getelementptr inbounds double, double* %7, i64 %79
  %108 = load double, double* %107, align 8, !tbaa !4
  %109 = fmul fast double %108, %106
  %110 = getelementptr inbounds double, double* %7, i64 %83
  %111 = load double, double* %110, align 8, !tbaa !4
  %112 = load double, double* %87, align 8, !tbaa !4
  %113 = fmul fast double %112, %111
  %114 = fadd fast double %113, %109
  %115 = load double, double* %91, align 8, !tbaa !4
  %116 = getelementptr inbounds double, double* %8, i64 %79
  %117 = load double, double* %116, align 8, !tbaa !4
  %118 = fmul fast double %117, %115
  %119 = fadd fast double %114, %118
  %120 = getelementptr inbounds double, double* %8, i64 %83
  %121 = load double, double* %120, align 8, !tbaa !4
  %122 = load double, double* %99, align 8, !tbaa !4
  %123 = fmul fast double %122, %121
  %124 = fadd fast double %119, %123
  %125 = getelementptr inbounds double, double* %1, i64 %104
  store double %124, double* %125, align 8, !tbaa !4
  %126 = load double, double* %77, align 8, !tbaa !4
  %127 = getelementptr inbounds double, double* %9, i64 %79
  %128 = load double, double* %127, align 8, !tbaa !4
  %129 = fmul fast double %128, %126
  %130 = getelementptr inbounds double, double* %9, i64 %83
  %131 = load double, double* %130, align 8, !tbaa !4
  %132 = load double, double* %87, align 8, !tbaa !4
  %133 = fmul fast double %132, %131
  %134 = fadd fast double %133, %129
  %135 = load double, double* %91, align 8, !tbaa !4
  %136 = getelementptr inbounds double, double* %10, i64 %79
  %137 = load double, double* %136, align 8, !tbaa !4
  %138 = fmul fast double %137, %135
  %139 = fadd fast double %134, %138
  %140 = getelementptr inbounds double, double* %10, i64 %83
  %141 = load double, double* %140, align 8, !tbaa !4
  %142 = load double, double* %99, align 8, !tbaa !4
  %143 = fmul fast double %142, %141
  %144 = fadd fast double %139, %143
  %145 = getelementptr inbounds double, double* %2, i64 %104
  store double %144, double* %145, align 8, !tbaa !4
  %146 = icmp eq i64 %61, 0
  br i1 %146, label %54, label %147

147:                                              ; preds = %60
  %148 = zext i32 %68 to i64
  br label %149

149:                                              ; preds = %149, %147
  %150 = phi i64 [ 1, %147 ], [ %231, %149 ]
  %151 = add nuw nsw i64 %150, %148
  %152 = mul nsw i64 %151, %25
  %153 = add nsw i64 %152, %52
  %154 = add nsw i64 %152, %53
  %155 = trunc i64 %151 to i32
  %156 = mul i32 %24, %155
  %157 = add i32 %49, %156
  %158 = add i32 %51, %156
  %159 = getelementptr inbounds double, double* %3, i64 %153
  %160 = load double, double* %159, align 8, !tbaa !4
  %161 = sext i32 %158 to i64
  %162 = getelementptr inbounds double, double* %5, i64 %161
  %163 = load double, double* %162, align 8, !tbaa !4
  %164 = fmul fast double %163, %160
  %165 = sext i32 %157 to i64
  %166 = getelementptr inbounds double, double* %5, i64 %165
  %167 = load double, double* %166, align 8, !tbaa !4
  %168 = getelementptr inbounds double, double* %3, i64 %154
  %169 = load double, double* %168, align 8, !tbaa !4
  %170 = fmul fast double %169, %167
  %171 = fadd fast double %170, %164
  %172 = getelementptr inbounds double, double* %4, i64 %153
  %173 = load double, double* %172, align 8, !tbaa !4
  %174 = getelementptr inbounds double, double* %6, i64 %161
  %175 = load double, double* %174, align 8, !tbaa !4
  %176 = fmul fast double %175, %173
  %177 = fadd fast double %171, %176
  %178 = getelementptr inbounds double, double* %6, i64 %165
  %179 = load double, double* %178, align 8, !tbaa !4
  %180 = getelementptr inbounds double, double* %4, i64 %154
  %181 = load double, double* %180, align 8, !tbaa !4
  %182 = fmul fast double %181, %179
  %183 = fadd fast double %177, %182
  %184 = fmul fast double %183, 2.000000e+00
  %185 = load double, double* %105, align 8, !tbaa !4
  %186 = fadd fast double %184, %185
  store double %186, double* %105, align 8, !tbaa !4
  %187 = load double, double* %159, align 8, !tbaa !4
  %188 = getelementptr inbounds double, double* %7, i64 %161
  %189 = load double, double* %188, align 8, !tbaa !4
  %190 = fmul fast double %189, %187
  %191 = getelementptr inbounds double, double* %7, i64 %165
  %192 = load double, double* %191, align 8, !tbaa !4
  %193 = load double, double* %168, align 8, !tbaa !4
  %194 = fmul fast double %193, %192
  %195 = fadd fast double %194, %190
  %196 = load double, double* %172, align 8, !tbaa !4
  %197 = getelementptr inbounds double, double* %8, i64 %161
  %198 = load double, double* %197, align 8, !tbaa !4
  %199 = fmul fast double %198, %196
  %200 = fadd fast double %195, %199
  %201 = getelementptr inbounds double, double* %8, i64 %165
  %202 = load double, double* %201, align 8, !tbaa !4
  %203 = load double, double* %180, align 8, !tbaa !4
  %204 = fmul fast double %203, %202
  %205 = fadd fast double %200, %204
  %206 = fmul fast double %205, 2.000000e+00
  %207 = load double, double* %125, align 8, !tbaa !4
  %208 = fadd fast double %206, %207
  store double %208, double* %125, align 8, !tbaa !4
  %209 = load double, double* %159, align 8, !tbaa !4
  %210 = getelementptr inbounds double, double* %9, i64 %161
  %211 = load double, double* %210, align 8, !tbaa !4
  %212 = fmul fast double %211, %209
  %213 = getelementptr inbounds double, double* %9, i64 %165
  %214 = load double, double* %213, align 8, !tbaa !4
  %215 = load double, double* %168, align 8, !tbaa !4
  %216 = fmul fast double %215, %214
  %217 = fadd fast double %216, %212
  %218 = load double, double* %172, align 8, !tbaa !4
  %219 = getelementptr inbounds double, double* %10, i64 %161
  %220 = load double, double* %219, align 8, !tbaa !4
  %221 = fmul fast double %220, %218
  %222 = fadd fast double %217, %221
  %223 = getelementptr inbounds double, double* %10, i64 %165
  %224 = load double, double* %223, align 8, !tbaa !4
  %225 = load double, double* %180, align 8, !tbaa !4
  %226 = fmul fast double %225, %224
  %227 = fadd fast double %222, %226
  %228 = fmul fast double %227, 2.000000e+00
  %229 = load double, double* %145, align 8, !tbaa !4
  %230 = fadd fast double %228, %229
  store double %230, double* %145, align 8, !tbaa !4
  %231 = add nuw nsw i64 %150, 1
  %232 = icmp eq i64 %231, %62
  br i1 %232, label %54, label %149
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z37cpuRadialSphericalHarmonicsPowerDerivIfEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_Piiii(float* %0, float* %1, float* %2, float* %3, float* %4, float* %5, float* %6, float* %7, float* %8, float* %9, float* %10, i32* %11, i32 %12, i32 %13, i32 %14) local_unnamed_addr #0 {
  %16 = add nsw i32 %13, 1
  %17 = mul nsw i32 %16, %13
  %18 = sdiv i32 %17, 2
  %19 = icmp sgt i32 %14, 0
  br i1 %19, label %20, label %35

20:                                               ; preds = %15
  %21 = icmp sgt i32 %17, 1
  %22 = icmp slt i32 %12, 0
  %23 = add i32 %12, 1
  %24 = mul i32 %14, %13
  %25 = sext i32 %13 to i64
  %26 = zext i32 %14 to i64
  %27 = sext i32 %23 to i64
  %28 = sext i32 %18 to i64
  %29 = zext i32 %18 to i64
  %30 = zext i32 %23 to i64
  br label %31

31:                                               ; preds = %36, %20
  %32 = phi i64 [ 0, %20 ], [ %37, %36 ]
  br i1 %21, label %33, label %36

33:                                               ; preds = %31
  %34 = trunc i64 %32 to i32
  br label %39

35:                                               ; preds = %36, %15
  ret void

36:                                               ; preds = %57, %31
  %37 = add nuw nsw i64 %32, 1
  %38 = icmp eq i64 %37, %26
  br i1 %38, label %35, label %31

39:                                               ; preds = %57, %33
  %40 = phi i64 [ 0, %33 ], [ %58, %57 ]
  %41 = getelementptr inbounds i32, i32* %11, i64 %40
  %42 = load i32, i32* %41, align 4, !tbaa !72
  %43 = add nsw i64 %40, %28
  %44 = getelementptr inbounds i32, i32* %11, i64 %43
  %45 = load i32, i32* %44, align 4, !tbaa !72
  br i1 %22, label %57, label %46

46:                                               ; preds = %39
  %47 = mul nsw i64 %40, %27
  %48 = mul nsw i32 %45, %14
  %49 = add i32 %48, %34
  %50 = mul nsw i32 %42, %14
  %51 = add i32 %50, %34
  %52 = sext i32 %45 to i64
  %53 = sext i32 %42 to i64
  br label %60

54:                                               ; preds = %149, %60
  %55 = add nuw nsw i64 %62, 1
  %56 = icmp eq i64 %64, %30
  br i1 %56, label %57, label %60

57:                                               ; preds = %54, %39
  %58 = add nuw nsw i64 %40, 1
  %59 = icmp eq i64 %58, %29
  br i1 %59, label %36, label %39

60:                                               ; preds = %54, %46
  %61 = phi i64 [ %64, %54 ], [ 0, %46 ]
  %62 = phi i64 [ %55, %54 ], [ 1, %46 ]
  %63 = phi i32 [ %65, %54 ], [ 0, %46 ]
  %64 = add nuw nsw i64 %61, 1
  %65 = add nuw nsw i32 %63, 1
  %66 = trunc i64 %61 to i32
  %67 = mul nsw i32 %65, %66
  %68 = lshr i32 %67, 1
  %69 = mul nsw i32 %68, %13
  %70 = add nsw i32 %69, %45
  %71 = add nsw i32 %69, %42
  %72 = add nsw i64 %61, %47
  %73 = mul i32 %24, %68
  %74 = add i32 %49, %73
  %75 = add i32 %51, %73
  %76 = sext i32 %70 to i64
  %77 = getelementptr inbounds float, float* %3, i64 %76
  %78 = load float, float* %77, align 4, !tbaa !8
  %79 = sext i32 %75 to i64
  %80 = getelementptr inbounds float, float* %5, i64 %79
  %81 = load float, float* %80, align 4, !tbaa !8
  %82 = fmul fast float %81, %78
  %83 = sext i32 %74 to i64
  %84 = getelementptr inbounds float, float* %5, i64 %83
  %85 = load float, float* %84, align 4, !tbaa !8
  %86 = sext i32 %71 to i64
  %87 = getelementptr inbounds float, float* %3, i64 %86
  %88 = load float, float* %87, align 4, !tbaa !8
  %89 = fmul fast float %88, %85
  %90 = fadd fast float %89, %82
  %91 = getelementptr inbounds float, float* %4, i64 %76
  %92 = load float, float* %91, align 4, !tbaa !8
  %93 = getelementptr inbounds float, float* %6, i64 %79
  %94 = load float, float* %93, align 4, !tbaa !8
  %95 = fmul fast float %94, %92
  %96 = fadd fast float %90, %95
  %97 = getelementptr inbounds float, float* %6, i64 %83
  %98 = load float, float* %97, align 4, !tbaa !8
  %99 = getelementptr inbounds float, float* %4, i64 %86
  %100 = load float, float* %99, align 4, !tbaa !8
  %101 = fmul fast float %100, %98
  %102 = fadd fast float %96, %101
  %103 = mul nsw i64 %72, %26
  %104 = add nsw i64 %103, %32
  %105 = getelementptr inbounds float, float* %0, i64 %104
  store float %102, float* %105, align 4, !tbaa !8
  %106 = load float, float* %77, align 4, !tbaa !8
  %107 = getelementptr inbounds float, float* %7, i64 %79
  %108 = load float, float* %107, align 4, !tbaa !8
  %109 = fmul fast float %108, %106
  %110 = getelementptr inbounds float, float* %7, i64 %83
  %111 = load float, float* %110, align 4, !tbaa !8
  %112 = load float, float* %87, align 4, !tbaa !8
  %113 = fmul fast float %112, %111
  %114 = fadd fast float %113, %109
  %115 = load float, float* %91, align 4, !tbaa !8
  %116 = getelementptr inbounds float, float* %8, i64 %79
  %117 = load float, float* %116, align 4, !tbaa !8
  %118 = fmul fast float %117, %115
  %119 = fadd fast float %114, %118
  %120 = getelementptr inbounds float, float* %8, i64 %83
  %121 = load float, float* %120, align 4, !tbaa !8
  %122 = load float, float* %99, align 4, !tbaa !8
  %123 = fmul fast float %122, %121
  %124 = fadd fast float %119, %123
  %125 = getelementptr inbounds float, float* %1, i64 %104
  store float %124, float* %125, align 4, !tbaa !8
  %126 = load float, float* %77, align 4, !tbaa !8
  %127 = getelementptr inbounds float, float* %9, i64 %79
  %128 = load float, float* %127, align 4, !tbaa !8
  %129 = fmul fast float %128, %126
  %130 = getelementptr inbounds float, float* %9, i64 %83
  %131 = load float, float* %130, align 4, !tbaa !8
  %132 = load float, float* %87, align 4, !tbaa !8
  %133 = fmul fast float %132, %131
  %134 = fadd fast float %133, %129
  %135 = load float, float* %91, align 4, !tbaa !8
  %136 = getelementptr inbounds float, float* %10, i64 %79
  %137 = load float, float* %136, align 4, !tbaa !8
  %138 = fmul fast float %137, %135
  %139 = fadd fast float %134, %138
  %140 = getelementptr inbounds float, float* %10, i64 %83
  %141 = load float, float* %140, align 4, !tbaa !8
  %142 = load float, float* %99, align 4, !tbaa !8
  %143 = fmul fast float %142, %141
  %144 = fadd fast float %139, %143
  %145 = getelementptr inbounds float, float* %2, i64 %104
  store float %144, float* %145, align 4, !tbaa !8
  %146 = icmp eq i64 %61, 0
  br i1 %146, label %54, label %147

147:                                              ; preds = %60
  %148 = zext i32 %68 to i64
  br label %149

149:                                              ; preds = %149, %147
  %150 = phi i64 [ 1, %147 ], [ %231, %149 ]
  %151 = add nuw nsw i64 %150, %148
  %152 = mul nsw i64 %151, %25
  %153 = add nsw i64 %152, %52
  %154 = add nsw i64 %152, %53
  %155 = trunc i64 %151 to i32
  %156 = mul i32 %24, %155
  %157 = add i32 %49, %156
  %158 = add i32 %51, %156
  %159 = getelementptr inbounds float, float* %3, i64 %153
  %160 = load float, float* %159, align 4, !tbaa !8
  %161 = sext i32 %158 to i64
  %162 = getelementptr inbounds float, float* %5, i64 %161
  %163 = load float, float* %162, align 4, !tbaa !8
  %164 = fmul fast float %163, %160
  %165 = sext i32 %157 to i64
  %166 = getelementptr inbounds float, float* %5, i64 %165
  %167 = load float, float* %166, align 4, !tbaa !8
  %168 = getelementptr inbounds float, float* %3, i64 %154
  %169 = load float, float* %168, align 4, !tbaa !8
  %170 = fmul fast float %169, %167
  %171 = fadd fast float %170, %164
  %172 = getelementptr inbounds float, float* %4, i64 %153
  %173 = load float, float* %172, align 4, !tbaa !8
  %174 = getelementptr inbounds float, float* %6, i64 %161
  %175 = load float, float* %174, align 4, !tbaa !8
  %176 = fmul fast float %175, %173
  %177 = fadd fast float %171, %176
  %178 = getelementptr inbounds float, float* %6, i64 %165
  %179 = load float, float* %178, align 4, !tbaa !8
  %180 = getelementptr inbounds float, float* %4, i64 %154
  %181 = load float, float* %180, align 4, !tbaa !8
  %182 = fmul fast float %181, %179
  %183 = fadd fast float %177, %182
  %184 = fmul fast float %183, 2.000000e+00
  %185 = load float, float* %105, align 4, !tbaa !8
  %186 = fadd fast float %184, %185
  store float %186, float* %105, align 4, !tbaa !8
  %187 = load float, float* %159, align 4, !tbaa !8
  %188 = getelementptr inbounds float, float* %7, i64 %161
  %189 = load float, float* %188, align 4, !tbaa !8
  %190 = fmul fast float %189, %187
  %191 = getelementptr inbounds float, float* %7, i64 %165
  %192 = load float, float* %191, align 4, !tbaa !8
  %193 = load float, float* %168, align 4, !tbaa !8
  %194 = fmul fast float %193, %192
  %195 = fadd fast float %194, %190
  %196 = load float, float* %172, align 4, !tbaa !8
  %197 = getelementptr inbounds float, float* %8, i64 %161
  %198 = load float, float* %197, align 4, !tbaa !8
  %199 = fmul fast float %198, %196
  %200 = fadd fast float %195, %199
  %201 = getelementptr inbounds float, float* %8, i64 %165
  %202 = load float, float* %201, align 4, !tbaa !8
  %203 = load float, float* %180, align 4, !tbaa !8
  %204 = fmul fast float %203, %202
  %205 = fadd fast float %200, %204
  %206 = fmul fast float %205, 2.000000e+00
  %207 = load float, float* %125, align 4, !tbaa !8
  %208 = fadd fast float %206, %207
  store float %208, float* %125, align 4, !tbaa !8
  %209 = load float, float* %159, align 4, !tbaa !8
  %210 = getelementptr inbounds float, float* %9, i64 %161
  %211 = load float, float* %210, align 4, !tbaa !8
  %212 = fmul fast float %211, %209
  %213 = getelementptr inbounds float, float* %9, i64 %165
  %214 = load float, float* %213, align 4, !tbaa !8
  %215 = load float, float* %168, align 4, !tbaa !8
  %216 = fmul fast float %215, %214
  %217 = fadd fast float %216, %212
  %218 = load float, float* %172, align 4, !tbaa !8
  %219 = getelementptr inbounds float, float* %10, i64 %161
  %220 = load float, float* %219, align 4, !tbaa !8
  %221 = fmul fast float %220, %218
  %222 = fadd fast float %217, %221
  %223 = getelementptr inbounds float, float* %10, i64 %165
  %224 = load float, float* %223, align 4, !tbaa !8
  %225 = load float, float* %180, align 4, !tbaa !8
  %226 = fmul fast float %225, %224
  %227 = fadd fast float %222, %226
  %228 = fmul fast float %227, 2.000000e+00
  %229 = load float, float* %145, align 4, !tbaa !8
  %230 = fadd fast float %228, %229
  store float %230, float* %145, align 4, !tbaa !8
  %231 = add nuw nsw i64 %150, 1
  %232 = icmp eq i64 %231, %62
  br i1 %232, label %54, label %149
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z37cpuRadialSphericalHarmonicsBispectrumIdEvPT_S1_S1_S1_PiS2_S2_S2_iii(double* %0, double* %1, double* %2, double* %3, i32* %4, i32* %5, i32* %6, i32* %7, i32 %8, i32 %9, i32 %10) local_unnamed_addr #0 {
  %12 = add nsw i32 %8, 1
  %13 = mul nsw i32 %12, %8
  %14 = sdiv i32 %13, 2
  %15 = icmp sgt i32 %13, 1
  br i1 %15, label %16, label %25

16:                                               ; preds = %11
  %17 = icmp sgt i32 %9, 0
  %18 = shl i32 %9, 1
  %19 = shl i32 %10, 1
  %20 = sext i32 %18 to i64
  %21 = sext i32 %9 to i64
  %22 = sext i32 %14 to i64
  %23 = zext i32 %14 to i64
  %24 = zext i32 %9 to i64
  br label %26

25:                                               ; preds = %36, %11
  ret void

26:                                               ; preds = %36, %16
  %27 = phi i64 [ 0, %16 ], [ %37, %36 ]
  %28 = getelementptr inbounds i32, i32* %4, i64 %27
  %29 = load i32, i32* %28, align 4, !tbaa !72
  %30 = add nsw i64 %27, %22
  %31 = getelementptr inbounds i32, i32* %4, i64 %30
  %32 = load i32, i32* %31, align 4, !tbaa !72
  br i1 %17, label %33, label %36

33:                                               ; preds = %26
  %34 = mul nsw i64 %27, %21
  %35 = load i32, i32* %7, align 4, !tbaa !72
  br label %39

36:                                               ; preds = %67, %26
  %37 = add nuw nsw i64 %27, 1
  %38 = icmp eq i64 %37, %23
  br i1 %38, label %25, label %26

39:                                               ; preds = %67, %33
  %40 = phi i32 [ %35, %33 ], [ %52, %67 ]
  %41 = phi i64 [ 0, %33 ], [ %50, %67 ]
  %42 = getelementptr inbounds i32, i32* %5, i64 %41
  %43 = load i32, i32* %42, align 4, !tbaa !72
  %44 = add nsw i64 %41, %21
  %45 = getelementptr inbounds i32, i32* %5, i64 %44
  %46 = load i32, i32* %45, align 4, !tbaa !72
  %47 = add nsw i64 %41, %20
  %48 = getelementptr inbounds i32, i32* %5, i64 %47
  %49 = load i32, i32* %48, align 4, !tbaa !72
  %50 = add nuw nsw i64 %41, 1
  %51 = getelementptr inbounds i32, i32* %7, i64 %50
  %52 = load i32, i32* %51, align 4, !tbaa !72
  %53 = sub i32 %52, %40
  %54 = icmp sgt i32 %53, 0
  br i1 %54, label %55, label %67

55:                                               ; preds = %39
  %56 = add nsw i32 %49, 1
  %57 = mul nsw i32 %56, %49
  %58 = sdiv i32 %57, 2
  %59 = add nsw i32 %46, 1
  %60 = mul nsw i32 %59, %46
  %61 = sdiv i32 %60, 2
  %62 = add nsw i32 %43, 1
  %63 = mul nsw i32 %62, %43
  %64 = sdiv i32 %63, 2
  %65 = sext i32 %40 to i64
  %66 = zext i32 %53 to i64
  br label %72

67:                                               ; preds = %166, %39
  %68 = phi double [ 0.000000e+00, %39 ], [ %181, %166 ]
  %69 = add nsw i64 %41, %34
  %70 = getelementptr inbounds double, double* %0, i64 %69
  store double %68, double* %70, align 8, !tbaa !4
  %71 = icmp eq i64 %50, %24
  br i1 %71, label %36, label %39

72:                                               ; preds = %166, %55
  %73 = phi i64 [ 0, %55 ], [ %182, %166 ]
  %74 = phi double [ 0.000000e+00, %55 ], [ %181, %166 ]
  %75 = add nsw i64 %73, %65
  %76 = getelementptr inbounds i32, i32* %6, i64 %75
  %77 = load i32, i32* %76, align 4, !tbaa !72
  %78 = trunc i64 %73 to i32
  %79 = add i32 %40, %78
  %80 = add i32 %79, %10
  %81 = sext i32 %80 to i64
  %82 = getelementptr inbounds i32, i32* %6, i64 %81
  %83 = load i32, i32* %82, align 4, !tbaa !72
  %84 = add i32 %79, %19
  %85 = sext i32 %84 to i64
  %86 = getelementptr inbounds i32, i32* %6, i64 %85
  %87 = load i32, i32* %86, align 4, !tbaa !72
  %88 = icmp sgt i32 %87, -1
  %89 = sub nsw i32 0, %87
  %90 = select i1 %88, i32 %87, i32 %89
  %91 = icmp sgt i32 %83, -1
  %92 = sub nsw i32 0, %83
  %93 = select i1 %91, i32 %83, i32 %92
  %94 = icmp sgt i32 %77, -1
  %95 = sub nsw i32 0, %77
  %96 = select i1 %94, i32 %77, i32 %95
  %97 = add nsw i32 %90, %58
  %98 = mul nsw i32 %97, %8
  %99 = add nsw i32 %98, %32
  %100 = add nsw i32 %93, %61
  %101 = mul nsw i32 %100, %8
  %102 = add nsw i32 %101, %29
  %103 = add nsw i32 %96, %64
  %104 = mul nsw i32 %103, %8
  %105 = add nsw i32 %104, %29
  %106 = and i32 %90, 1
  %107 = icmp eq i32 %106, 0
  %108 = select fast i1 %107, double 1.000000e+00, double -1.000000e+00
  %109 = and i32 %93, 1
  %110 = icmp eq i32 %109, 0
  %111 = select fast i1 %110, double 1.000000e+00, double -1.000000e+00
  %112 = and i32 %96, 1
  %113 = icmp eq i32 %112, 0
  %114 = select fast i1 %113, double 1.000000e+00, double -1.000000e+00
  br i1 %88, label %115, label %121

115:                                              ; preds = %72
  %116 = sext i32 %99 to i64
  %117 = getelementptr inbounds double, double* %1, i64 %116
  %118 = load double, double* %117, align 8, !tbaa !4
  %119 = getelementptr inbounds double, double* %2, i64 %116
  %120 = load double, double* %119, align 8, !tbaa !4
  br label %130

121:                                              ; preds = %72
  %122 = fsub fast double -0.000000e+00, %108
  %123 = sext i32 %99 to i64
  %124 = getelementptr inbounds double, double* %1, i64 %123
  %125 = load double, double* %124, align 8, !tbaa !4
  %126 = fmul fast double %125, %122
  %127 = getelementptr inbounds double, double* %2, i64 %123
  %128 = load double, double* %127, align 8, !tbaa !4
  %129 = fmul fast double %128, %108
  br label %130

130:                                              ; preds = %121, %115
  %131 = phi double [ %118, %115 ], [ %126, %121 ]
  %132 = phi double [ %120, %115 ], [ %129, %121 ]
  br i1 %91, label %133, label %139

133:                                              ; preds = %130
  %134 = sext i32 %102 to i64
  %135 = getelementptr inbounds double, double* %1, i64 %134
  %136 = load double, double* %135, align 8, !tbaa !4
  %137 = getelementptr inbounds double, double* %2, i64 %134
  %138 = load double, double* %137, align 8, !tbaa !4
  br label %148

139:                                              ; preds = %130
  %140 = fsub fast double -0.000000e+00, %111
  %141 = sext i32 %102 to i64
  %142 = getelementptr inbounds double, double* %1, i64 %141
  %143 = load double, double* %142, align 8, !tbaa !4
  %144 = fmul fast double %143, %140
  %145 = getelementptr inbounds double, double* %2, i64 %141
  %146 = load double, double* %145, align 8, !tbaa !4
  %147 = fmul fast double %146, %111
  br label %148

148:                                              ; preds = %139, %133
  %149 = phi double [ %136, %133 ], [ %144, %139 ]
  %150 = phi double [ %138, %133 ], [ %147, %139 ]
  br i1 %94, label %151, label %157

151:                                              ; preds = %148
  %152 = sext i32 %105 to i64
  %153 = getelementptr inbounds double, double* %1, i64 %152
  %154 = load double, double* %153, align 8, !tbaa !4
  %155 = getelementptr inbounds double, double* %2, i64 %152
  %156 = load double, double* %155, align 8, !tbaa !4
  br label %166

157:                                              ; preds = %148
  %158 = fsub fast double -0.000000e+00, %114
  %159 = sext i32 %105 to i64
  %160 = getelementptr inbounds double, double* %1, i64 %159
  %161 = load double, double* %160, align 8, !tbaa !4
  %162 = fmul fast double %161, %158
  %163 = getelementptr inbounds double, double* %2, i64 %159
  %164 = load double, double* %163, align 8, !tbaa !4
  %165 = fmul fast double %164, %114
  br label %166

166:                                              ; preds = %157, %151
  %167 = phi double [ %154, %151 ], [ %162, %157 ]
  %168 = phi double [ %156, %151 ], [ %165, %157 ]
  %169 = getelementptr inbounds double, double* %3, i64 %75
  %170 = load double, double* %169, align 8, !tbaa !4
  %171 = fmul fast double %167, %149
  %172 = fmul fast double %168, %149
  %173 = fmul fast double %167, %150
  %174 = fmul fast double %168, %150
  %175 = fadd fast double %173, %172
  %176 = fmul fast double %175, %132
  %177 = fsub fast double %171, %174
  %178 = fmul fast double %177, %131
  %179 = fadd fast double %178, %176
  %180 = fmul fast double %179, %170
  %181 = fadd fast double %180, %74
  %182 = add nuw nsw i64 %73, 1
  %183 = icmp eq i64 %182, %66
  br i1 %183, label %67, label %72
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z37cpuRadialSphericalHarmonicsBispectrumIfEvPT_S1_S1_S1_PiS2_S2_S2_iii(float* %0, float* %1, float* %2, float* %3, i32* %4, i32* %5, i32* %6, i32* %7, i32 %8, i32 %9, i32 %10) local_unnamed_addr #0 {
  %12 = add nsw i32 %8, 1
  %13 = mul nsw i32 %12, %8
  %14 = sdiv i32 %13, 2
  %15 = icmp sgt i32 %13, 1
  br i1 %15, label %16, label %25

16:                                               ; preds = %11
  %17 = icmp sgt i32 %9, 0
  %18 = shl i32 %9, 1
  %19 = shl i32 %10, 1
  %20 = sext i32 %18 to i64
  %21 = sext i32 %9 to i64
  %22 = sext i32 %14 to i64
  %23 = zext i32 %14 to i64
  %24 = zext i32 %9 to i64
  br label %26

25:                                               ; preds = %36, %11
  ret void

26:                                               ; preds = %36, %16
  %27 = phi i64 [ 0, %16 ], [ %37, %36 ]
  %28 = getelementptr inbounds i32, i32* %4, i64 %27
  %29 = load i32, i32* %28, align 4, !tbaa !72
  %30 = add nsw i64 %27, %22
  %31 = getelementptr inbounds i32, i32* %4, i64 %30
  %32 = load i32, i32* %31, align 4, !tbaa !72
  br i1 %17, label %33, label %36

33:                                               ; preds = %26
  %34 = mul nsw i64 %27, %21
  %35 = load i32, i32* %7, align 4, !tbaa !72
  br label %39

36:                                               ; preds = %67, %26
  %37 = add nuw nsw i64 %27, 1
  %38 = icmp eq i64 %37, %23
  br i1 %38, label %25, label %26

39:                                               ; preds = %67, %33
  %40 = phi i32 [ %35, %33 ], [ %52, %67 ]
  %41 = phi i64 [ 0, %33 ], [ %50, %67 ]
  %42 = getelementptr inbounds i32, i32* %5, i64 %41
  %43 = load i32, i32* %42, align 4, !tbaa !72
  %44 = add nsw i64 %41, %21
  %45 = getelementptr inbounds i32, i32* %5, i64 %44
  %46 = load i32, i32* %45, align 4, !tbaa !72
  %47 = add nsw i64 %41, %20
  %48 = getelementptr inbounds i32, i32* %5, i64 %47
  %49 = load i32, i32* %48, align 4, !tbaa !72
  %50 = add nuw nsw i64 %41, 1
  %51 = getelementptr inbounds i32, i32* %7, i64 %50
  %52 = load i32, i32* %51, align 4, !tbaa !72
  %53 = sub i32 %52, %40
  %54 = icmp sgt i32 %53, 0
  br i1 %54, label %55, label %67

55:                                               ; preds = %39
  %56 = add nsw i32 %49, 1
  %57 = mul nsw i32 %56, %49
  %58 = sdiv i32 %57, 2
  %59 = add nsw i32 %46, 1
  %60 = mul nsw i32 %59, %46
  %61 = sdiv i32 %60, 2
  %62 = add nsw i32 %43, 1
  %63 = mul nsw i32 %62, %43
  %64 = sdiv i32 %63, 2
  %65 = sext i32 %40 to i64
  %66 = zext i32 %53 to i64
  br label %72

67:                                               ; preds = %166, %39
  %68 = phi float [ 0.000000e+00, %39 ], [ %181, %166 ]
  %69 = add nsw i64 %41, %34
  %70 = getelementptr inbounds float, float* %0, i64 %69
  store float %68, float* %70, align 4, !tbaa !8
  %71 = icmp eq i64 %50, %24
  br i1 %71, label %36, label %39

72:                                               ; preds = %166, %55
  %73 = phi i64 [ 0, %55 ], [ %182, %166 ]
  %74 = phi float [ 0.000000e+00, %55 ], [ %181, %166 ]
  %75 = add nsw i64 %73, %65
  %76 = getelementptr inbounds i32, i32* %6, i64 %75
  %77 = load i32, i32* %76, align 4, !tbaa !72
  %78 = trunc i64 %73 to i32
  %79 = add i32 %40, %78
  %80 = add i32 %79, %10
  %81 = sext i32 %80 to i64
  %82 = getelementptr inbounds i32, i32* %6, i64 %81
  %83 = load i32, i32* %82, align 4, !tbaa !72
  %84 = add i32 %79, %19
  %85 = sext i32 %84 to i64
  %86 = getelementptr inbounds i32, i32* %6, i64 %85
  %87 = load i32, i32* %86, align 4, !tbaa !72
  %88 = icmp sgt i32 %87, -1
  %89 = sub nsw i32 0, %87
  %90 = select i1 %88, i32 %87, i32 %89
  %91 = icmp sgt i32 %83, -1
  %92 = sub nsw i32 0, %83
  %93 = select i1 %91, i32 %83, i32 %92
  %94 = icmp sgt i32 %77, -1
  %95 = sub nsw i32 0, %77
  %96 = select i1 %94, i32 %77, i32 %95
  %97 = add nsw i32 %90, %58
  %98 = mul nsw i32 %97, %8
  %99 = add nsw i32 %98, %32
  %100 = add nsw i32 %93, %61
  %101 = mul nsw i32 %100, %8
  %102 = add nsw i32 %101, %29
  %103 = add nsw i32 %96, %64
  %104 = mul nsw i32 %103, %8
  %105 = add nsw i32 %104, %29
  %106 = and i32 %90, 1
  %107 = icmp eq i32 %106, 0
  %108 = select i1 %107, float 1.000000e+00, float -1.000000e+00
  %109 = and i32 %93, 1
  %110 = icmp eq i32 %109, 0
  %111 = select i1 %110, float 1.000000e+00, float -1.000000e+00
  %112 = and i32 %96, 1
  %113 = icmp eq i32 %112, 0
  %114 = select i1 %113, float 1.000000e+00, float -1.000000e+00
  br i1 %88, label %115, label %121

115:                                              ; preds = %72
  %116 = sext i32 %99 to i64
  %117 = getelementptr inbounds float, float* %1, i64 %116
  %118 = load float, float* %117, align 4, !tbaa !8
  %119 = getelementptr inbounds float, float* %2, i64 %116
  %120 = load float, float* %119, align 4, !tbaa !8
  br label %130

121:                                              ; preds = %72
  %122 = fsub fast float -0.000000e+00, %108
  %123 = sext i32 %99 to i64
  %124 = getelementptr inbounds float, float* %1, i64 %123
  %125 = load float, float* %124, align 4, !tbaa !8
  %126 = fmul fast float %125, %122
  %127 = getelementptr inbounds float, float* %2, i64 %123
  %128 = load float, float* %127, align 4, !tbaa !8
  %129 = fmul fast float %128, %108
  br label %130

130:                                              ; preds = %121, %115
  %131 = phi float [ %118, %115 ], [ %126, %121 ]
  %132 = phi float [ %120, %115 ], [ %129, %121 ]
  br i1 %91, label %133, label %139

133:                                              ; preds = %130
  %134 = sext i32 %102 to i64
  %135 = getelementptr inbounds float, float* %1, i64 %134
  %136 = load float, float* %135, align 4, !tbaa !8
  %137 = getelementptr inbounds float, float* %2, i64 %134
  %138 = load float, float* %137, align 4, !tbaa !8
  br label %148

139:                                              ; preds = %130
  %140 = fsub fast float -0.000000e+00, %111
  %141 = sext i32 %102 to i64
  %142 = getelementptr inbounds float, float* %1, i64 %141
  %143 = load float, float* %142, align 4, !tbaa !8
  %144 = fmul fast float %143, %140
  %145 = getelementptr inbounds float, float* %2, i64 %141
  %146 = load float, float* %145, align 4, !tbaa !8
  %147 = fmul fast float %146, %111
  br label %148

148:                                              ; preds = %139, %133
  %149 = phi float [ %136, %133 ], [ %144, %139 ]
  %150 = phi float [ %138, %133 ], [ %147, %139 ]
  br i1 %94, label %151, label %157

151:                                              ; preds = %148
  %152 = sext i32 %105 to i64
  %153 = getelementptr inbounds float, float* %1, i64 %152
  %154 = load float, float* %153, align 4, !tbaa !8
  %155 = getelementptr inbounds float, float* %2, i64 %152
  %156 = load float, float* %155, align 4, !tbaa !8
  br label %166

157:                                              ; preds = %148
  %158 = fsub fast float -0.000000e+00, %114
  %159 = sext i32 %105 to i64
  %160 = getelementptr inbounds float, float* %1, i64 %159
  %161 = load float, float* %160, align 4, !tbaa !8
  %162 = fmul fast float %161, %158
  %163 = getelementptr inbounds float, float* %2, i64 %159
  %164 = load float, float* %163, align 4, !tbaa !8
  %165 = fmul fast float %164, %114
  br label %166

166:                                              ; preds = %157, %151
  %167 = phi float [ %154, %151 ], [ %162, %157 ]
  %168 = phi float [ %156, %151 ], [ %165, %157 ]
  %169 = getelementptr inbounds float, float* %3, i64 %75
  %170 = load float, float* %169, align 4, !tbaa !8
  %171 = fmul fast float %167, %149
  %172 = fmul fast float %168, %149
  %173 = fmul fast float %167, %150
  %174 = fmul fast float %168, %150
  %175 = fadd fast float %173, %172
  %176 = fmul fast float %175, %132
  %177 = fsub fast float %171, %174
  %178 = fmul fast float %177, %131
  %179 = fadd fast float %178, %176
  %180 = fmul fast float %179, %170
  %181 = fadd fast float %180, %74
  %182 = add nuw nsw i64 %73, 1
  %183 = icmp eq i64 %182, %66
  br i1 %183, label %67, label %72
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z42cpuRadialSphericalHarmonicsBispectrumDerivIdEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_PiS2_S2_S2_iiii(double* %0, double* %1, double* %2, double* %3, double* %4, double* %5, double* %6, double* %7, double* %8, double* %9, double* %10, double* %11, i32* %12, i32* %13, i32* %14, i32* %15, i32 %16, i32 %17, i32 %18, i32 %19) local_unnamed_addr #0 {
  %21 = add nsw i32 %17, 1
  %22 = mul nsw i32 %21, %17
  %23 = sdiv i32 %22, 2
  %24 = icmp sgt i32 %19, 0
  br i1 %24, label %25, label %41

25:                                               ; preds = %20
  %26 = icmp sgt i32 %22, 1
  %27 = icmp sgt i32 %16, 0
  %28 = shl i32 %16, 1
  %29 = shl i32 %18, 1
  %30 = mul i32 %19, %17
  %31 = zext i32 %19 to i64
  %32 = sext i32 %28 to i64
  %33 = sext i32 %16 to i64
  %34 = sext i32 %23 to i64
  %35 = zext i32 %23 to i64
  %36 = zext i32 %16 to i64
  br label %37

37:                                               ; preds = %42, %25
  %38 = phi i64 [ 0, %25 ], [ %43, %42 ]
  br i1 %26, label %39, label %42

39:                                               ; preds = %37
  %40 = trunc i64 %38 to i32
  br label %45

41:                                               ; preds = %42, %20
  ret void

42:                                               ; preds = %58, %37
  %43 = add nuw nsw i64 %38, 1
  %44 = icmp eq i64 %43, %31
  br i1 %44, label %41, label %37

45:                                               ; preds = %58, %39
  %46 = phi i64 [ 0, %39 ], [ %59, %58 ]
  %47 = getelementptr inbounds i32, i32* %12, i64 %46
  %48 = load i32, i32* %47, align 4, !tbaa !72
  %49 = add nsw i64 %46, %34
  %50 = getelementptr inbounds i32, i32* %12, i64 %49
  %51 = load i32, i32* %50, align 4, !tbaa !72
  br i1 %27, label %52, label %58

52:                                               ; preds = %45
  %53 = mul nsw i64 %46, %33
  %54 = mul nsw i32 %48, %19
  %55 = add i32 %54, %40
  br label %61

56:                                               ; preds = %303, %61
  %57 = icmp eq i64 %77, %36
  br i1 %57, label %58, label %61

58:                                               ; preds = %56, %45
  %59 = add nuw nsw i64 %46, 1
  %60 = icmp eq i64 %59, %35
  br i1 %60, label %42, label %45

61:                                               ; preds = %56, %52
  %62 = phi i64 [ 0, %52 ], [ %77, %56 ]
  %63 = getelementptr inbounds i32, i32* %13, i64 %62
  %64 = load i32, i32* %63, align 4, !tbaa !72
  %65 = add nsw i64 %62, %33
  %66 = getelementptr inbounds i32, i32* %13, i64 %65
  %67 = load i32, i32* %66, align 4, !tbaa !72
  %68 = add nsw i64 %62, %32
  %69 = getelementptr inbounds i32, i32* %13, i64 %68
  %70 = load i32, i32* %69, align 4, !tbaa !72
  %71 = add nsw i64 %62, %53
  %72 = mul nsw i64 %71, %31
  %73 = add nsw i64 %72, %38
  %74 = getelementptr inbounds double, double* %0, i64 %73
  store double 0.000000e+00, double* %74, align 8, !tbaa !4
  %75 = getelementptr inbounds double, double* %1, i64 %73
  store double 0.000000e+00, double* %75, align 8, !tbaa !4
  %76 = getelementptr inbounds double, double* %2, i64 %73
  store double 0.000000e+00, double* %76, align 8, !tbaa !4
  %77 = add nuw nsw i64 %62, 1
  %78 = getelementptr inbounds i32, i32* %15, i64 %77
  %79 = load i32, i32* %78, align 4, !tbaa !72
  %80 = getelementptr inbounds i32, i32* %15, i64 %62
  %81 = load i32, i32* %80, align 4, !tbaa !72
  %82 = sub i32 %79, %81
  %83 = icmp sgt i32 %82, 0
  br i1 %83, label %84, label %56

84:                                               ; preds = %61
  %85 = add nsw i32 %70, 1
  %86 = mul nsw i32 %85, %70
  %87 = sdiv i32 %86, 2
  %88 = add nsw i32 %67, 1
  %89 = mul nsw i32 %88, %67
  %90 = sdiv i32 %89, 2
  %91 = add nsw i32 %64, 1
  %92 = mul nsw i32 %91, %64
  %93 = sdiv i32 %92, 2
  %94 = sext i32 %81 to i64
  %95 = zext i32 %82 to i64
  br label %96

96:                                               ; preds = %303, %84
  %97 = phi i64 [ 0, %84 ], [ %400, %303 ]
  %98 = add nsw i64 %97, %94
  %99 = getelementptr inbounds i32, i32* %14, i64 %98
  %100 = load i32, i32* %99, align 4, !tbaa !72
  %101 = trunc i64 %97 to i32
  %102 = add i32 %81, %101
  %103 = add i32 %102, %18
  %104 = sext i32 %103 to i64
  %105 = getelementptr inbounds i32, i32* %14, i64 %104
  %106 = load i32, i32* %105, align 4, !tbaa !72
  %107 = add i32 %102, %29
  %108 = sext i32 %107 to i64
  %109 = getelementptr inbounds i32, i32* %14, i64 %108
  %110 = load i32, i32* %109, align 4, !tbaa !72
  %111 = icmp sgt i32 %110, -1
  %112 = sub nsw i32 0, %110
  %113 = select i1 %111, i32 %110, i32 %112
  %114 = icmp sgt i32 %106, -1
  %115 = sub nsw i32 0, %106
  %116 = select i1 %114, i32 %106, i32 %115
  %117 = icmp sgt i32 %100, -1
  %118 = sub nsw i32 0, %100
  %119 = select i1 %117, i32 %100, i32 %118
  %120 = add nsw i32 %113, %87
  %121 = mul i32 %120, %17
  %122 = add i32 %121, %51
  %123 = add nsw i32 %116, %90
  %124 = mul nsw i32 %123, %17
  %125 = add nsw i32 %124, %48
  %126 = add nsw i32 %119, %93
  %127 = mul nsw i32 %126, %17
  %128 = add nsw i32 %127, %48
  %129 = mul i32 %122, %19
  %130 = add i32 %129, %40
  %131 = mul i32 %30, %123
  %132 = add i32 %55, %131
  %133 = mul i32 %30, %126
  %134 = add i32 %55, %133
  %135 = and i32 %113, 1
  %136 = icmp eq i32 %135, 0
  %137 = select fast i1 %136, double 1.000000e+00, double -1.000000e+00
  %138 = and i32 %116, 1
  %139 = icmp eq i32 %138, 0
  %140 = select fast i1 %139, double 1.000000e+00, double -1.000000e+00
  %141 = and i32 %119, 1
  %142 = icmp eq i32 %141, 0
  %143 = select fast i1 %142, double 1.000000e+00, double -1.000000e+00
  br i1 %111, label %144, label %163

144:                                              ; preds = %96
  %145 = sext i32 %122 to i64
  %146 = getelementptr inbounds double, double* %3, i64 %145
  %147 = load double, double* %146, align 8, !tbaa !4
  %148 = getelementptr inbounds double, double* %4, i64 %145
  %149 = load double, double* %148, align 8, !tbaa !4
  %150 = sext i32 %130 to i64
  %151 = getelementptr inbounds double, double* %5, i64 %150
  %152 = load double, double* %151, align 8, !tbaa !4
  %153 = getelementptr inbounds double, double* %7, i64 %150
  %154 = load double, double* %153, align 8, !tbaa !4
  %155 = getelementptr inbounds double, double* %9, i64 %150
  %156 = load double, double* %155, align 8, !tbaa !4
  %157 = getelementptr inbounds double, double* %6, i64 %150
  %158 = load double, double* %157, align 8, !tbaa !4
  %159 = getelementptr inbounds double, double* %8, i64 %150
  %160 = load double, double* %159, align 8, !tbaa !4
  %161 = getelementptr inbounds double, double* %10, i64 %150
  %162 = load double, double* %161, align 8, !tbaa !4
  br label %191

163:                                              ; preds = %96
  %164 = fneg fast double %137
  %165 = sext i32 %122 to i64
  %166 = getelementptr inbounds double, double* %3, i64 %165
  %167 = load double, double* %166, align 8, !tbaa !4
  %168 = fmul fast double %167, %164
  %169 = getelementptr inbounds double, double* %4, i64 %165
  %170 = load double, double* %169, align 8, !tbaa !4
  %171 = fmul fast double %170, %137
  %172 = sext i32 %130 to i64
  %173 = getelementptr inbounds double, double* %5, i64 %172
  %174 = load double, double* %173, align 8, !tbaa !4
  %175 = fmul fast double %174, %164
  %176 = getelementptr inbounds double, double* %7, i64 %172
  %177 = load double, double* %176, align 8, !tbaa !4
  %178 = fmul fast double %177, %164
  %179 = getelementptr inbounds double, double* %9, i64 %172
  %180 = load double, double* %179, align 8, !tbaa !4
  %181 = fmul fast double %180, %164
  %182 = getelementptr inbounds double, double* %6, i64 %172
  %183 = load double, double* %182, align 8, !tbaa !4
  %184 = fmul fast double %183, %137
  %185 = getelementptr inbounds double, double* %8, i64 %172
  %186 = load double, double* %185, align 8, !tbaa !4
  %187 = fmul fast double %186, %137
  %188 = getelementptr inbounds double, double* %10, i64 %172
  %189 = load double, double* %188, align 8, !tbaa !4
  %190 = fmul fast double %189, %137
  br label %191

191:                                              ; preds = %163, %144
  %192 = phi double [ %149, %144 ], [ %171, %163 ]
  %193 = phi double [ %147, %144 ], [ %168, %163 ]
  %194 = phi double [ %152, %144 ], [ %175, %163 ]
  %195 = phi double [ %158, %144 ], [ %184, %163 ]
  %196 = phi double [ %154, %144 ], [ %178, %163 ]
  %197 = phi double [ %160, %144 ], [ %187, %163 ]
  %198 = phi double [ %156, %144 ], [ %181, %163 ]
  %199 = phi double [ %162, %144 ], [ %190, %163 ]
  br i1 %114, label %200, label %219

200:                                              ; preds = %191
  %201 = sext i32 %125 to i64
  %202 = getelementptr inbounds double, double* %3, i64 %201
  %203 = load double, double* %202, align 8, !tbaa !4
  %204 = getelementptr inbounds double, double* %4, i64 %201
  %205 = load double, double* %204, align 8, !tbaa !4
  %206 = sext i32 %132 to i64
  %207 = getelementptr inbounds double, double* %5, i64 %206
  %208 = load double, double* %207, align 8, !tbaa !4
  %209 = getelementptr inbounds double, double* %7, i64 %206
  %210 = load double, double* %209, align 8, !tbaa !4
  %211 = getelementptr inbounds double, double* %9, i64 %206
  %212 = load double, double* %211, align 8, !tbaa !4
  %213 = getelementptr inbounds double, double* %6, i64 %206
  %214 = load double, double* %213, align 8, !tbaa !4
  %215 = getelementptr inbounds double, double* %8, i64 %206
  %216 = load double, double* %215, align 8, !tbaa !4
  %217 = getelementptr inbounds double, double* %10, i64 %206
  %218 = load double, double* %217, align 8, !tbaa !4
  br label %247

219:                                              ; preds = %191
  %220 = fneg fast double %140
  %221 = sext i32 %125 to i64
  %222 = getelementptr inbounds double, double* %3, i64 %221
  %223 = load double, double* %222, align 8, !tbaa !4
  %224 = fmul fast double %223, %220
  %225 = getelementptr inbounds double, double* %4, i64 %221
  %226 = load double, double* %225, align 8, !tbaa !4
  %227 = fmul fast double %226, %140
  %228 = sext i32 %132 to i64
  %229 = getelementptr inbounds double, double* %5, i64 %228
  %230 = load double, double* %229, align 8, !tbaa !4
  %231 = fmul fast double %230, %220
  %232 = getelementptr inbounds double, double* %7, i64 %228
  %233 = load double, double* %232, align 8, !tbaa !4
  %234 = fmul fast double %233, %220
  %235 = getelementptr inbounds double, double* %9, i64 %228
  %236 = load double, double* %235, align 8, !tbaa !4
  %237 = fmul fast double %236, %220
  %238 = getelementptr inbounds double, double* %6, i64 %228
  %239 = load double, double* %238, align 8, !tbaa !4
  %240 = fmul fast double %239, %140
  %241 = getelementptr inbounds double, double* %8, i64 %228
  %242 = load double, double* %241, align 8, !tbaa !4
  %243 = fmul fast double %242, %140
  %244 = getelementptr inbounds double, double* %10, i64 %228
  %245 = load double, double* %244, align 8, !tbaa !4
  %246 = fmul fast double %245, %140
  br label %247

247:                                              ; preds = %219, %200
  %248 = phi double [ %203, %200 ], [ %224, %219 ]
  %249 = phi double [ %205, %200 ], [ %227, %219 ]
  %250 = phi double [ %208, %200 ], [ %231, %219 ]
  %251 = phi double [ %214, %200 ], [ %240, %219 ]
  %252 = phi double [ %210, %200 ], [ %234, %219 ]
  %253 = phi double [ %216, %200 ], [ %243, %219 ]
  %254 = phi double [ %212, %200 ], [ %237, %219 ]
  %255 = phi double [ %218, %200 ], [ %246, %219 ]
  br i1 %117, label %256, label %275

256:                                              ; preds = %247
  %257 = sext i32 %128 to i64
  %258 = getelementptr inbounds double, double* %3, i64 %257
  %259 = load double, double* %258, align 8, !tbaa !4
  %260 = getelementptr inbounds double, double* %4, i64 %257
  %261 = load double, double* %260, align 8, !tbaa !4
  %262 = sext i32 %134 to i64
  %263 = getelementptr inbounds double, double* %5, i64 %262
  %264 = load double, double* %263, align 8, !tbaa !4
  %265 = getelementptr inbounds double, double* %7, i64 %262
  %266 = load double, double* %265, align 8, !tbaa !4
  %267 = getelementptr inbounds double, double* %9, i64 %262
  %268 = load double, double* %267, align 8, !tbaa !4
  %269 = getelementptr inbounds double, double* %6, i64 %262
  %270 = load double, double* %269, align 8, !tbaa !4
  %271 = getelementptr inbounds double, double* %8, i64 %262
  %272 = load double, double* %271, align 8, !tbaa !4
  %273 = getelementptr inbounds double, double* %10, i64 %262
  %274 = load double, double* %273, align 8, !tbaa !4
  br label %303

275:                                              ; preds = %247
  %276 = fneg fast double %143
  %277 = sext i32 %128 to i64
  %278 = getelementptr inbounds double, double* %3, i64 %277
  %279 = load double, double* %278, align 8, !tbaa !4
  %280 = fmul fast double %279, %276
  %281 = getelementptr inbounds double, double* %4, i64 %277
  %282 = load double, double* %281, align 8, !tbaa !4
  %283 = fmul fast double %282, %143
  %284 = sext i32 %134 to i64
  %285 = getelementptr inbounds double, double* %5, i64 %284
  %286 = load double, double* %285, align 8, !tbaa !4
  %287 = fmul fast double %286, %276
  %288 = getelementptr inbounds double, double* %7, i64 %284
  %289 = load double, double* %288, align 8, !tbaa !4
  %290 = fmul fast double %289, %276
  %291 = getelementptr inbounds double, double* %9, i64 %284
  %292 = load double, double* %291, align 8, !tbaa !4
  %293 = fmul fast double %292, %276
  %294 = getelementptr inbounds double, double* %6, i64 %284
  %295 = load double, double* %294, align 8, !tbaa !4
  %296 = fmul fast double %295, %143
  %297 = getelementptr inbounds double, double* %8, i64 %284
  %298 = load double, double* %297, align 8, !tbaa !4
  %299 = fmul fast double %298, %143
  %300 = getelementptr inbounds double, double* %10, i64 %284
  %301 = load double, double* %300, align 8, !tbaa !4
  %302 = fmul fast double %301, %143
  br label %303

303:                                              ; preds = %275, %256
  %304 = phi double [ %259, %256 ], [ %280, %275 ]
  %305 = phi double [ %261, %256 ], [ %283, %275 ]
  %306 = phi double [ %264, %256 ], [ %287, %275 ]
  %307 = phi double [ %270, %256 ], [ %296, %275 ]
  %308 = phi double [ %266, %256 ], [ %290, %275 ]
  %309 = phi double [ %272, %256 ], [ %299, %275 ]
  %310 = phi double [ %268, %256 ], [ %293, %275 ]
  %311 = phi double [ %274, %256 ], [ %302, %275 ]
  %312 = getelementptr inbounds double, double* %11, i64 %98
  %313 = load double, double* %312, align 8, !tbaa !4
  %314 = fmul fast double %248, %194
  %315 = fmul fast double %250, %193
  %316 = fadd fast double %315, %314
  %317 = fmul fast double %304, %316
  %318 = fmul fast double %248, %193
  %319 = fmul fast double %306, %318
  %320 = fmul fast double %250, %192
  %321 = fmul fast double %248, %195
  %322 = fadd fast double %320, %321
  %323 = fmul fast double %248, %192
  %324 = fmul fast double %306, %192
  %325 = fmul fast double %304, %195
  %326 = fadd fast double %324, %325
  %327 = fmul fast double %326, %249
  %328 = fmul fast double %304, %192
  %329 = fmul fast double %328, %251
  %330 = fmul fast double %249, %194
  %331 = fmul fast double %251, %193
  %332 = fadd fast double %331, %330
  %333 = fmul fast double %249, %193
  %334 = load double, double* %74, align 8, !tbaa !4
  %335 = fsub fast double %323, %333
  %336 = fmul fast double %335, %307
  %337 = fsub fast double %322, %332
  %338 = fmul fast double %337, %305
  %339 = fadd fast double %329, %317
  %340 = fadd fast double %339, %319
  %341 = fadd fast double %340, %327
  %342 = fadd fast double %341, %338
  %343 = fadd fast double %342, %336
  %344 = fmul fast double %313, %343
  %345 = fadd fast double %334, %344
  store double %345, double* %74, align 8, !tbaa !4
  %346 = fmul fast double %248, %196
  %347 = fmul fast double %252, %193
  %348 = fadd fast double %347, %346
  %349 = fmul fast double %304, %348
  %350 = fmul fast double %308, %318
  %351 = fmul fast double %252, %192
  %352 = fmul fast double %248, %197
  %353 = fadd fast double %351, %352
  %354 = fmul fast double %308, %192
  %355 = fmul fast double %304, %197
  %356 = fadd fast double %354, %355
  %357 = fmul fast double %356, %249
  %358 = fmul fast double %328, %253
  %359 = fmul fast double %249, %196
  %360 = fmul fast double %253, %193
  %361 = fadd fast double %360, %359
  %362 = load double, double* %75, align 8, !tbaa !4
  %363 = fmul fast double %335, %309
  %364 = fsub fast double %353, %361
  %365 = fmul fast double %364, %305
  %366 = fadd fast double %358, %349
  %367 = fadd fast double %366, %350
  %368 = fadd fast double %367, %365
  %369 = fadd fast double %368, %357
  %370 = fadd fast double %369, %363
  %371 = fmul fast double %370, %313
  %372 = fadd fast double %362, %371
  store double %372, double* %75, align 8, !tbaa !4
  %373 = fmul fast double %248, %198
  %374 = fmul fast double %254, %193
  %375 = fadd fast double %374, %373
  %376 = fmul fast double %304, %375
  %377 = fmul fast double %310, %318
  %378 = fmul fast double %254, %192
  %379 = fmul fast double %248, %199
  %380 = fadd fast double %378, %379
  %381 = fmul fast double %310, %192
  %382 = fmul fast double %304, %199
  %383 = fadd fast double %381, %382
  %384 = fmul fast double %383, %249
  %385 = fmul fast double %328, %255
  %386 = fmul fast double %249, %198
  %387 = fmul fast double %255, %193
  %388 = fadd fast double %387, %386
  %389 = load double, double* %76, align 8, !tbaa !4
  %390 = fmul fast double %335, %311
  %391 = fsub fast double %380, %388
  %392 = fmul fast double %391, %305
  %393 = fadd fast double %385, %376
  %394 = fadd fast double %393, %392
  %395 = fadd fast double %394, %377
  %396 = fadd fast double %395, %384
  %397 = fadd fast double %396, %390
  %398 = fmul fast double %397, %313
  %399 = fadd fast double %389, %398
  store double %399, double* %76, align 8, !tbaa !4
  %400 = add nuw nsw i64 %97, 1
  %401 = icmp eq i64 %400, %95
  br i1 %401, label %56, label %96
}
; Function Attrs: nounwind ssp uwtable
define weak_odr void @_Z42cpuRadialSphericalHarmonicsBispectrumDerivIfEvPT_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_S1_PiS2_S2_S2_iiii(float* %0, float* %1, float* %2, float* %3, float* %4, float* %5, float* %6, float* %7, float* %8, float* %9, float* %10, float* %11, i32* %12, i32* %13, i32* %14, i32* %15, i32 %16, i32 %17, i32 %18, i32 %19) local_unnamed_addr #0 {
  %21 = add nsw i32 %17, 1
  %22 = mul nsw i32 %21, %17
  %23 = sdiv i32 %22, 2
  %24 = icmp sgt i32 %19, 0
  br i1 %24, label %25, label %41

25:                                               ; preds = %20
  %26 = icmp sgt i32 %22, 1
  %27 = icmp sgt i32 %16, 0
  %28 = shl i32 %16, 1
  %29 = shl i32 %18, 1
  %30 = mul i32 %19, %17
  %31 = zext i32 %19 to i64
  %32 = sext i32 %28 to i64
  %33 = sext i32 %16 to i64
  %34 = sext i32 %23 to i64
  %35 = zext i32 %23 to i64
  %36 = zext i32 %16 to i64
  br label %37

37:                                               ; preds = %42, %25
  %38 = phi i64 [ 0, %25 ], [ %43, %42 ]
  br i1 %26, label %39, label %42

39:                                               ; preds = %37
  %40 = trunc i64 %38 to i32
  br label %45

41:                                               ; preds = %42, %20
  ret void

42:                                               ; preds = %58, %37
  %43 = add nuw nsw i64 %38, 1
  %44 = icmp eq i64 %43, %31
  br i1 %44, label %41, label %37

45:                                               ; preds = %58, %39
  %46 = phi i64 [ 0, %39 ], [ %59, %58 ]
  %47 = getelementptr inbounds i32, i32* %12, i64 %46
  %48 = load i32, i32* %47, align 4, !tbaa !72
  %49 = add nsw i64 %46, %34
  %50 = getelementptr inbounds i32, i32* %12, i64 %49
  %51 = load i32, i32* %50, align 4, !tbaa !72
  br i1 %27, label %52, label %58

52:                                               ; preds = %45
  %53 = mul nsw i64 %46, %33
  %54 = mul nsw i32 %48, %19
  %55 = add i32 %54, %40
  br label %61

56:                                               ; preds = %303, %61
  %57 = icmp eq i64 %77, %36
  br i1 %57, label %58, label %61

58:                                               ; preds = %56, %45
  %59 = add nuw nsw i64 %46, 1
  %60 = icmp eq i64 %59, %35
  br i1 %60, label %42, label %45

61:                                               ; preds = %56, %52
  %62 = phi i64 [ 0, %52 ], [ %77, %56 ]
  %63 = getelementptr inbounds i32, i32* %13, i64 %62
  %64 = load i32, i32* %63, align 4, !tbaa !72
  %65 = add nsw i64 %62, %33
  %66 = getelementptr inbounds i32, i32* %13, i64 %65
  %67 = load i32, i32* %66, align 4, !tbaa !72
  %68 = add nsw i64 %62, %32
  %69 = getelementptr inbounds i32, i32* %13, i64 %68
  %70 = load i32, i32* %69, align 4, !tbaa !72
  %71 = add nsw i64 %62, %53
  %72 = mul nsw i64 %71, %31
  %73 = add nsw i64 %72, %38
  %74 = getelementptr inbounds float, float* %0, i64 %73
  store float 0.000000e+00, float* %74, align 4, !tbaa !8
  %75 = getelementptr inbounds float, float* %1, i64 %73
  store float 0.000000e+00, float* %75, align 4, !tbaa !8
  %76 = getelementptr inbounds float, float* %2, i64 %73
  store float 0.000000e+00, float* %76, align 4, !tbaa !8
  %77 = add nuw nsw i64 %62, 1
  %78 = getelementptr inbounds i32, i32* %15, i64 %77
  %79 = load i32, i32* %78, align 4, !tbaa !72
  %80 = getelementptr inbounds i32, i32* %15, i64 %62
  %81 = load i32, i32* %80, align 4, !tbaa !72
  %82 = sub i32 %79, %81
  %83 = icmp sgt i32 %82, 0
  br i1 %83, label %84, label %56

84:                                               ; preds = %61
  %85 = add nsw i32 %70, 1
  %86 = mul nsw i32 %85, %70
  %87 = sdiv i32 %86, 2
  %88 = add nsw i32 %67, 1
  %89 = mul nsw i32 %88, %67
  %90 = sdiv i32 %89, 2
  %91 = add nsw i32 %64, 1
  %92 = mul nsw i32 %91, %64
  %93 = sdiv i32 %92, 2
  %94 = sext i32 %81 to i64
  %95 = zext i32 %82 to i64
  br label %96

96:                                               ; preds = %303, %84
  %97 = phi i64 [ 0, %84 ], [ %400, %303 ]
  %98 = add nsw i64 %97, %94
  %99 = getelementptr inbounds i32, i32* %14, i64 %98
  %100 = load i32, i32* %99, align 4, !tbaa !72
  %101 = trunc i64 %97 to i32
  %102 = add i32 %81, %101
  %103 = add i32 %102, %18
  %104 = sext i32 %103 to i64
  %105 = getelementptr inbounds i32, i32* %14, i64 %104
  %106 = load i32, i32* %105, align 4, !tbaa !72
  %107 = add i32 %102, %29
  %108 = sext i32 %107 to i64
  %109 = getelementptr inbounds i32, i32* %14, i64 %108
  %110 = load i32, i32* %109, align 4, !tbaa !72
  %111 = icmp sgt i32 %110, -1
  %112 = sub nsw i32 0, %110
  %113 = select i1 %111, i32 %110, i32 %112
  %114 = icmp sgt i32 %106, -1
  %115 = sub nsw i32 0, %106
  %116 = select i1 %114, i32 %106, i32 %115
  %117 = icmp sgt i32 %100, -1
  %118 = sub nsw i32 0, %100
  %119 = select i1 %117, i32 %100, i32 %118
  %120 = add nsw i32 %113, %87
  %121 = mul i32 %120, %17
  %122 = add i32 %121, %51
  %123 = add nsw i32 %116, %90
  %124 = mul nsw i32 %123, %17
  %125 = add nsw i32 %124, %48
  %126 = add nsw i32 %119, %93
  %127 = mul nsw i32 %126, %17
  %128 = add nsw i32 %127, %48
  %129 = mul i32 %122, %19
  %130 = add i32 %129, %40
  %131 = mul i32 %30, %123
  %132 = add i32 %55, %131
  %133 = mul i32 %30, %126
  %134 = add i32 %55, %133
  %135 = and i32 %113, 1
  %136 = icmp eq i32 %135, 0
  %137 = select i1 %136, float 1.000000e+00, float -1.000000e+00
  %138 = and i32 %116, 1
  %139 = icmp eq i32 %138, 0
  %140 = select i1 %139, float 1.000000e+00, float -1.000000e+00
  %141 = and i32 %119, 1
  %142 = icmp eq i32 %141, 0
  %143 = select i1 %142, float 1.000000e+00, float -1.000000e+00
  br i1 %111, label %144, label %163

144:                                              ; preds = %96
  %145 = sext i32 %122 to i64
  %146 = getelementptr inbounds float, float* %3, i64 %145
  %147 = load float, float* %146, align 4, !tbaa !8
  %148 = getelementptr inbounds float, float* %4, i64 %145
  %149 = load float, float* %148, align 4, !tbaa !8
  %150 = sext i32 %130 to i64
  %151 = getelementptr inbounds float, float* %5, i64 %150
  %152 = load float, float* %151, align 4, !tbaa !8
  %153 = getelementptr inbounds float, float* %7, i64 %150
  %154 = load float, float* %153, align 4, !tbaa !8
  %155 = getelementptr inbounds float, float* %9, i64 %150
  %156 = load float, float* %155, align 4, !tbaa !8
  %157 = getelementptr inbounds float, float* %6, i64 %150
  %158 = load float, float* %157, align 4, !tbaa !8
  %159 = getelementptr inbounds float, float* %8, i64 %150
  %160 = load float, float* %159, align 4, !tbaa !8
  %161 = getelementptr inbounds float, float* %10, i64 %150
  %162 = load float, float* %161, align 4, !tbaa !8
  br label %191

163:                                              ; preds = %96
  %164 = fneg fast float %137
  %165 = sext i32 %122 to i64
  %166 = getelementptr inbounds float, float* %3, i64 %165
  %167 = load float, float* %166, align 4, !tbaa !8
  %168 = fmul fast float %167, %164
  %169 = getelementptr inbounds float, float* %4, i64 %165
  %170 = load float, float* %169, align 4, !tbaa !8
  %171 = fmul fast float %170, %137
  %172 = sext i32 %130 to i64
  %173 = getelementptr inbounds float, float* %5, i64 %172
  %174 = load float, float* %173, align 4, !tbaa !8
  %175 = fmul fast float %174, %164
  %176 = getelementptr inbounds float, float* %7, i64 %172
  %177 = load float, float* %176, align 4, !tbaa !8
  %178 = fmul fast float %177, %164
  %179 = getelementptr inbounds float, float* %9, i64 %172
  %180 = load float, float* %179, align 4, !tbaa !8
  %181 = fmul fast float %180, %164
  %182 = getelementptr inbounds float, float* %6, i64 %172
  %183 = load float, float* %182, align 4, !tbaa !8
  %184 = fmul fast float %183, %137
  %185 = getelementptr inbounds float, float* %8, i64 %172
  %186 = load float, float* %185, align 4, !tbaa !8
  %187 = fmul fast float %186, %137
  %188 = getelementptr inbounds float, float* %10, i64 %172
  %189 = load float, float* %188, align 4, !tbaa !8
  %190 = fmul fast float %189, %137
  br label %191

191:                                              ; preds = %163, %144
  %192 = phi float [ %149, %144 ], [ %171, %163 ]
  %193 = phi float [ %147, %144 ], [ %168, %163 ]
  %194 = phi float [ %152, %144 ], [ %175, %163 ]
  %195 = phi float [ %158, %144 ], [ %184, %163 ]
  %196 = phi float [ %154, %144 ], [ %178, %163 ]
  %197 = phi float [ %160, %144 ], [ %187, %163 ]
  %198 = phi float [ %156, %144 ], [ %181, %163 ]
  %199 = phi float [ %162, %144 ], [ %190, %163 ]
  br i1 %114, label %200, label %219

200:                                              ; preds = %191
  %201 = sext i32 %125 to i64
  %202 = getelementptr inbounds float, float* %3, i64 %201
  %203 = load float, float* %202, align 4, !tbaa !8
  %204 = getelementptr inbounds float, float* %4, i64 %201
  %205 = load float, float* %204, align 4, !tbaa !8
  %206 = sext i32 %132 to i64
  %207 = getelementptr inbounds float, float* %5, i64 %206
  %208 = load float, float* %207, align 4, !tbaa !8
  %209 = getelementptr inbounds float, float* %7, i64 %206
  %210 = load float, float* %209, align 4, !tbaa !8
  %211 = getelementptr inbounds float, float* %9, i64 %206
  %212 = load float, float* %211, align 4, !tbaa !8
  %213 = getelementptr inbounds float, float* %6, i64 %206
  %214 = load float, float* %213, align 4, !tbaa !8
  %215 = getelementptr inbounds float, float* %8, i64 %206
  %216 = load float, float* %215, align 4, !tbaa !8
  %217 = getelementptr inbounds float, float* %10, i64 %206
  %218 = load float, float* %217, align 4, !tbaa !8
  br label %247

219:                                              ; preds = %191
  %220 = fneg fast float %140
  %221 = sext i32 %125 to i64
  %222 = getelementptr inbounds float, float* %3, i64 %221
  %223 = load float, float* %222, align 4, !tbaa !8
  %224 = fmul fast float %223, %220
  %225 = getelementptr inbounds float, float* %4, i64 %221
  %226 = load float, float* %225, align 4, !tbaa !8
  %227 = fmul fast float %226, %140
  %228 = sext i32 %132 to i64
  %229 = getelementptr inbounds float, float* %5, i64 %228
  %230 = load float, float* %229, align 4, !tbaa !8
  %231 = fmul fast float %230, %220
  %232 = getelementptr inbounds float, float* %7, i64 %228
  %233 = load float, float* %232, align 4, !tbaa !8
  %234 = fmul fast float %233, %220
  %235 = getelementptr inbounds float, float* %9, i64 %228
  %236 = load float, float* %235, align 4, !tbaa !8
  %237 = fmul fast float %236, %220
  %238 = getelementptr inbounds float, float* %6, i64 %228
  %239 = load float, float* %238, align 4, !tbaa !8
  %240 = fmul fast float %239, %140
  %241 = getelementptr inbounds float, float* %8, i64 %228
  %242 = load float, float* %241, align 4, !tbaa !8
  %243 = fmul fast float %242, %140
  %244 = getelementptr inbounds float, float* %10, i64 %228
  %245 = load float, float* %244, align 4, !tbaa !8
  %246 = fmul fast float %245, %140
  br label %247

247:                                              ; preds = %219, %200
  %248 = phi float [ %203, %200 ], [ %224, %219 ]
  %249 = phi float [ %205, %200 ], [ %227, %219 ]
  %250 = phi float [ %208, %200 ], [ %231, %219 ]
  %251 = phi float [ %214, %200 ], [ %240, %219 ]
  %252 = phi float [ %210, %200 ], [ %234, %219 ]
  %253 = phi float [ %216, %200 ], [ %243, %219 ]
  %254 = phi float [ %212, %200 ], [ %237, %219 ]
  %255 = phi float [ %218, %200 ], [ %246, %219 ]
  br i1 %117, label %256, label %275

256:                                              ; preds = %247
  %257 = sext i32 %128 to i64
  %258 = getelementptr inbounds float, float* %3, i64 %257
  %259 = load float, float* %258, align 4, !tbaa !8
  %260 = getelementptr inbounds float, float* %4, i64 %257
  %261 = load float, float* %260, align 4, !tbaa !8
  %262 = sext i32 %134 to i64
  %263 = getelementptr inbounds float, float* %5, i64 %262
  %264 = load float, float* %263, align 4, !tbaa !8
  %265 = getelementptr inbounds float, float* %7, i64 %262
  %266 = load float, float* %265, align 4, !tbaa !8
  %267 = getelementptr inbounds float, float* %9, i64 %262
  %268 = load float, float* %267, align 4, !tbaa !8
  %269 = getelementptr inbounds float, float* %6, i64 %262
  %270 = load float, float* %269, align 4, !tbaa !8
  %271 = getelementptr inbounds float, float* %8, i64 %262
  %272 = load float, float* %271, align 4, !tbaa !8
  %273 = getelementptr inbounds float, float* %10, i64 %262
  %274 = load float, float* %273, align 4, !tbaa !8
  br label %303

275:                                              ; preds = %247
  %276 = fneg fast float %143
  %277 = sext i32 %128 to i64
  %278 = getelementptr inbounds float, float* %3, i64 %277
  %279 = load float, float* %278, align 4, !tbaa !8
  %280 = fmul fast float %279, %276
  %281 = getelementptr inbounds float, float* %4, i64 %277
  %282 = load float, float* %281, align 4, !tbaa !8
  %283 = fmul fast float %282, %143
  %284 = sext i32 %134 to i64
  %285 = getelementptr inbounds float, float* %5, i64 %284
  %286 = load float, float* %285, align 4, !tbaa !8
  %287 = fmul fast float %286, %276
  %288 = getelementptr inbounds float, float* %7, i64 %284
  %289 = load float, float* %288, align 4, !tbaa !8
  %290 = fmul fast float %289, %276
  %291 = getelementptr inbounds float, float* %9, i64 %284
  %292 = load float, float* %291, align 4, !tbaa !8
  %293 = fmul fast float %292, %276
  %294 = getelementptr inbounds float, float* %6, i64 %284
  %295 = load float, float* %294, align 4, !tbaa !8
  %296 = fmul fast float %295, %143
  %297 = getelementptr inbounds float, float* %8, i64 %284
  %298 = load float, float* %297, align 4, !tbaa !8
  %299 = fmul fast float %298, %143
  %300 = getelementptr inbounds float, float* %10, i64 %284
  %301 = load float, float* %300, align 4, !tbaa !8
  %302 = fmul fast float %301, %143
  br label %303

303:                                              ; preds = %275, %256
  %304 = phi float [ %259, %256 ], [ %280, %275 ]
  %305 = phi float [ %261, %256 ], [ %283, %275 ]
  %306 = phi float [ %264, %256 ], [ %287, %275 ]
  %307 = phi float [ %270, %256 ], [ %296, %275 ]
  %308 = phi float [ %266, %256 ], [ %290, %275 ]
  %309 = phi float [ %272, %256 ], [ %299, %275 ]
  %310 = phi float [ %268, %256 ], [ %293, %275 ]
  %311 = phi float [ %274, %256 ], [ %302, %275 ]
  %312 = getelementptr inbounds float, float* %11, i64 %98
  %313 = load float, float* %312, align 4, !tbaa !8
  %314 = fmul fast float %248, %194
  %315 = fmul fast float %250, %193
  %316 = fadd fast float %315, %314
  %317 = fmul fast float %304, %316
  %318 = fmul fast float %248, %193
  %319 = fmul fast float %306, %318
  %320 = fmul fast float %250, %192
  %321 = fmul fast float %248, %195
  %322 = fadd fast float %320, %321
  %323 = fmul fast float %248, %192
  %324 = fmul fast float %306, %192
  %325 = fmul fast float %304, %195
  %326 = fadd fast float %324, %325
  %327 = fmul fast float %326, %249
  %328 = fmul fast float %304, %192
  %329 = fmul fast float %328, %251
  %330 = fmul fast float %249, %194
  %331 = fmul fast float %251, %193
  %332 = fadd fast float %331, %330
  %333 = fmul fast float %249, %193
  %334 = load float, float* %74, align 4, !tbaa !8
  %335 = fsub fast float %323, %333
  %336 = fmul fast float %335, %307
  %337 = fsub fast float %322, %332
  %338 = fmul fast float %337, %305
  %339 = fadd fast float %329, %317
  %340 = fadd fast float %339, %319
  %341 = fadd fast float %340, %327
  %342 = fadd fast float %341, %338
  %343 = fadd fast float %342, %336
  %344 = fmul fast float %313, %343
  %345 = fadd fast float %334, %344
  store float %345, float* %74, align 4, !tbaa !8
  %346 = fmul fast float %248, %196
  %347 = fmul fast float %252, %193
  %348 = fadd fast float %347, %346
  %349 = fmul fast float %304, %348
  %350 = fmul fast float %308, %318
  %351 = fmul fast float %252, %192
  %352 = fmul fast float %248, %197
  %353 = fadd fast float %351, %352
  %354 = fmul fast float %308, %192
  %355 = fmul fast float %304, %197
  %356 = fadd fast float %354, %355
  %357 = fmul fast float %356, %249
  %358 = fmul fast float %328, %253
  %359 = fmul fast float %249, %196
  %360 = fmul fast float %253, %193
  %361 = fadd fast float %360, %359
  %362 = load float, float* %75, align 4, !tbaa !8
  %363 = fmul fast float %335, %309
  %364 = fsub fast float %353, %361
  %365 = fmul fast float %364, %305
  %366 = fadd fast float %358, %349
  %367 = fadd fast float %366, %350
  %368 = fadd fast float %367, %365
  %369 = fadd fast float %368, %357
  %370 = fadd fast float %369, %363
  %371 = fmul fast float %370, %313
  %372 = fadd fast float %362, %371
  store float %372, float* %75, align 4, !tbaa !8
  %373 = fmul fast float %248, %198
  %374 = fmul fast float %254, %193
  %375 = fadd fast float %374, %373
  %376 = fmul fast float %304, %375
  %377 = fmul fast float %310, %318
  %378 = fmul fast float %254, %192
  %379 = fmul fast float %248, %199
  %380 = fadd fast float %378, %379
  %381 = fmul fast float %310, %192
  %382 = fmul fast float %304, %199
  %383 = fadd fast float %381, %382
  %384 = fmul fast float %383, %249
  %385 = fmul fast float %328, %255
  %386 = fmul fast float %249, %198
  %387 = fmul fast float %255, %193
  %388 = fadd fast float %387, %386
  %389 = load float, float* %76, align 4, !tbaa !8
  %390 = fmul fast float %335, %311
  %391 = fsub fast float %380, %388
  %392 = fmul fast float %391, %305
  %393 = fadd fast float %385, %376
  %394 = fadd fast float %393, %392
  %395 = fadd fast float %394, %377
  %396 = fadd fast float %395, %384
  %397 = fadd fast float %396, %390
  %398 = fmul fast float %397, %313
  %399 = fadd fast float %389, %398
  store float %399, float* %76, align 4, !tbaa !8
  %400 = add nuw nsw i64 %97, 1
  %401 = icmp eq i64 %400, %95
  br i1 %401, label %56, label %96
}
; Function Attrs: norecurse ssp uwtable
define i32 @main() local_unnamed_addr #3 {
  %1 = alloca [30 x double], align 16
  %2 = bitcast [30 x double]* %1 to i8*
  %3 = alloca [30 x double], align 16
  %4 = bitcast [30 x double]* %3 to i8*
  %5 = alloca [30 x double], align 16
  %6 = alloca [30 x double], align 16
  %7 = alloca [30 x double], align 16
  %8 = alloca [30 x double], align 16
  %9 = alloca [30 x double], align 16
  %10 = alloca [30 x double], align 16
  %11 = alloca [30 x double], align 16
  %12 = alloca [30 x double], align 16
  %13 = alloca [30 x double], align 16
  %14 = alloca [30 x double], align 16
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %2) #8
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %4) #8
  %15 = bitcast [30 x double]* %5 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %15) #8
  call void @llvm.memset.p0i8.i64(i8* nonnull align 16 dereferenceable(240) %15, i8 0, i64 240, i1 false)
  %16 = bitcast [30 x double]* %6 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %16) #8
  call void @llvm.memset.p0i8.i64(i8* nonnull align 16 dereferenceable(240) %16, i8 0, i64 240, i1 false)
  %17 = bitcast [30 x double]* %7 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %17) #8
  %18 = bitcast [30 x double]* %8 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %18) #8
  %19 = bitcast [30 x double]* %9 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %19) #8
  %20 = bitcast [30 x double]* %10 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %20) #8
  %21 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 0
  %22 = getelementptr inbounds [30 x double], [30 x double]* %5, i64 0, i64 0
  %23 = getelementptr inbounds [30 x double], [30 x double]* %7, i64 0, i64 0
  %24 = getelementptr inbounds [30 x double], [30 x double]* %9, i64 0, i64 0
  %25 = getelementptr inbounds [30 x double], [30 x double]* %10, i64 0, i64 0
  %26 = tail call fast double @llvm.sqrt.f64(double 0x7FF8000000000000) #8
  call void @llvm.memset.p0i8.i64(i8* nonnull align 16 dereferenceable(160) %4, i8 0, i64 160, i1 false)
  call void @memset_pattern16(i8* nonnull %2, i8* bitcast ([2 x double]* @.memset_pattern to i8*), i64 80)
  %27 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 10
  store double %26, double* %27, align 16, !tbaa !4
  %28 = getelementptr inbounds [30 x double], [30 x double]* %6, i64 0, i64 0
  %29 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 20
  %30 = getelementptr inbounds [30 x double], [30 x double]* %3, i64 0, i64 20
  %31 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 11
  store double %26, double* %31, align 8, !tbaa !4
  %32 = bitcast double* %29 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %32, align 16, !tbaa !4
  %33 = bitcast double* %30 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %33, align 16, !tbaa !4
  %34 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 12
  store double %26, double* %34, align 16, !tbaa !4
  %35 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 22
  %36 = getelementptr inbounds [30 x double], [30 x double]* %3, i64 0, i64 22
  %37 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 13
  store double %26, double* %37, align 8, !tbaa !4
  %38 = bitcast double* %35 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %38, align 16, !tbaa !4
  %39 = bitcast double* %36 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %39, align 16, !tbaa !4
  %40 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 14
  store double %26, double* %40, align 16, !tbaa !4
  %41 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 24
  %42 = getelementptr inbounds [30 x double], [30 x double]* %3, i64 0, i64 24
  %43 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 15
  store double %26, double* %43, align 8, !tbaa !4
  %44 = bitcast double* %41 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %44, align 16, !tbaa !4
  %45 = bitcast double* %42 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %45, align 16, !tbaa !4
  %46 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 16
  store double %26, double* %46, align 16, !tbaa !4
  %47 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 26
  %48 = getelementptr inbounds [30 x double], [30 x double]* %3, i64 0, i64 26
  %49 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 17
  store double %26, double* %49, align 8, !tbaa !4
  %50 = bitcast double* %47 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %50, align 16, !tbaa !4
  %51 = bitcast double* %48 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %51, align 16, !tbaa !4
  %52 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 18
  store double %26, double* %52, align 16, !tbaa !4
  %53 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 28
  %54 = getelementptr inbounds [30 x double], [30 x double]* %3, i64 0, i64 28
  %55 = getelementptr inbounds [30 x double], [30 x double]* %1, i64 0, i64 19
  store double %26, double* %55, align 8, !tbaa !4
  %56 = bitcast double* %53 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %56, align 16, !tbaa !4
  %57 = bitcast double* %54 to <2 x double>*
  store <2 x double> zeroinitializer, <2 x double>* %57, align 16, !tbaa !4
  %58 = getelementptr inbounds [30 x double], [30 x double]* %3, i64 0, i64 0
  %59 = getelementptr inbounds [30 x double], [30 x double]* %8, i64 0, i64 0
  %60 = getelementptr inbounds [30 x double], [30 x double]* %10, i64 0, i64 1
  %61 = bitcast [30 x double]* %7 to <2 x double>*
  store <2 x double> <double 1.000000e+00, double -0.000000e+00>, <2 x double>* %61, align 16, !tbaa !4
  store double %26, double* %25, align 16, !tbaa !4
  store double %26, double* %60, align 8, !tbaa !4
  %62 = load double, double* %21, align 16, !tbaa !4
  %63 = tail call i32 (i8*, ...) @printf(i8* nonnull dereferenceable(1) getelementptr inbounds ([12 x i8], [12 x i8]* @.str, i64 0, i64 0), double %62)
  %64 = bitcast [30 x double]* %11 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %64) #8
  call void @llvm.memset.p0i8.i64(i8* nonnull align 16 dereferenceable(240) %64, i8 0, i64 240, i1 false)
  %65 = bitcast [30 x double]* %12 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %65) #8
  call void @llvm.memset.p0i8.i64(i8* nonnull align 16 dereferenceable(240) %65, i8 0, i64 240, i1 false)
  %66 = bitcast [30 x double]* %13 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %66) #8
  call void @llvm.memset.p0i8.i64(i8* nonnull align 16 dereferenceable(240) %66, i8 0, i64 240, i1 false)
  %67 = getelementptr inbounds [30 x double], [30 x double]* %13, i64 0, i64 0
  store double 1.000000e+00, double* %67, align 16
  %68 = bitcast [30 x double]* %14 to i8*
  call void @llvm.lifetime.start.p0i8(i64 240, i8* nonnull %68) #8
  call void @llvm.memset.p0i8.i64(i8* nonnull align 16 dereferenceable(240) %68, i8 0, i64 240, i1 false)
  %69 = getelementptr inbounds [30 x double], [30 x double]* %11, i64 0, i64 0
  %70 = getelementptr inbounds [30 x double], [30 x double]* %12, i64 0, i64 0
  %71 = getelementptr inbounds [30 x double], [30 x double]* %14, i64 0, i64 0
  %72 = load i32, i32* @enzyme_const, align 4, !tbaa !72
  call void @_Z17__enzyme_autodiffIJPdS0_S0_S0_S0_S0_S0_S0_iS0_iS0_iS0_iS0_idiiiiEEvPvDpT_(i8* bitcast (void (double*, double*, double*, double*, double*, double*, double*, double*, double, i32, i32)* @_Z21cpuSphericalHarmonicsIdEvPT_S1_S1_S1_S1_S1_S1_S1_S0_ii to i8*), double* nonnull %21, double* nonnull %69, double* nonnull %58, double* nonnull %70, double* nonnull %22, double* nonnull %67, double* nonnull %28, double* nonnull %71, i32 %72, double* nonnull %23, i32 %72, double* nonnull %59, i32 %72, double* nonnull %24, i32 %72, double* nonnull %25, i32 %72, double 0x400921FB54442D18, i32 %72, i32 1, i32 %72, i32 10)
  %73 = load double, double* %69, align 16, !tbaa !4
  %74 = call i32 (i8*, ...) @printf(i8* nonnull dereferenceable(1) getelementptr inbounds ([14 x i8], [14 x i8]* @.str.1, i64 0, i64 0), double %73)
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %68) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %66) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %65) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %64) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %20) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %19) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %18) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %17) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %16) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %15) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %4) #8
  call void @llvm.lifetime.end.p0i8(i64 240, i8* nonnull %2) #8
  ret i32 0
}
; Function Attrs: argmemonly nounwind willreturn writeonly
declare void @llvm.memset.p0i8.i64(i8* nocapture writeonly, i8, i64, i1 immarg) #4
; Function Attrs: nofree nounwind
declare i32 @printf(i8* nocapture readonly, ...) local_unnamed_addr #5
declare void @_Z17__enzyme_autodiffIJPdS0_S0_S0_S0_S0_S0_S0_iS0_iS0_iS0_iS0_idiiiiEEvPvDpT_(i8*, double*, double*, double*, double*, double*, double*, double*, double*, i32, double*, i32, double*, i32, double*, i32, double*, i32, double, i32, i32, i32, i32) local_unnamed_addr #6
; Function Attrs: nounwind readnone speculatable willreturn
declare float @llvm.sqrt.f32(float) #2
; Function Attrs: nounwind readnone speculatable willreturn
declare float @llvm.cos.f32(float) #2
; Function Attrs: nounwind readnone speculatable willreturn
declare float @llvm.sin.f32(float) #2
; Function Attrs: argmemonly nofree
declare void @memset_pattern16(i8* nocapture, i8* nocapture readonly, i64) local_unnamed_addr #7

attributes #0 = { nounwind ssp uwtable "correctly-rounded-divide-sqrt-fp-math"="false" "darwin-stkchk-strong-link" "disable-tail-calls"="false" "frame-pointer"="all" "less-precise-fpmad"="false" "min-legal-vector-width"="0" "no-infs-fp-math"="true" "no-jump-tables"="false" "no-nans-fp-math"="true" "no-signed-zeros-fp-math"="true" "no-trapping-math"="true" "probe-stack"="___chkstk_darwin" "stack-protector-buffer-size"="8" "target-cpu"="penryn" "target-features"="+cx16,+cx8,+fxsr,+mmx,+sahf,+sse,+sse2,+sse3,+sse4.1,+ssse3,+x87" "unsafe-fp-math"="true" "use-soft-float"="false" }
attributes #1 = { argmemonly nounwind willreturn }
attributes #2 = { nounwind readnone speculatable willreturn }
attributes #3 = { norecurse ssp uwtable "correctly-rounded-divide-sqrt-fp-math"="false" "darwin-stkchk-strong-link" "disable-tail-calls"="false" "frame-pointer"="all" "less-precise-fpmad"="false" "min-legal-vector-width"="0" "no-infs-fp-math"="true" "no-jump-tables"="false" "no-nans-fp-math"="true" "no-signed-zeros-fp-math"="true" "no-trapping-math"="true" "probe-stack"="___chkstk_darwin" "stack-protector-buffer-size"="8" "target-cpu"="penryn" "target-features"="+cx16,+cx8,+fxsr,+mmx,+sahf,+sse,+sse2,+sse3,+sse4.1,+ssse3,+x87" "unsafe-fp-math"="true" "use-soft-float"="false" }
attributes #4 = { argmemonly nounwind willreturn writeonly }
attributes #5 = { nofree nounwind "correctly-rounded-divide-sqrt-fp-math"="false" "darwin-stkchk-strong-link" "disable-tail-calls"="false" "frame-pointer"="all" "less-precise-fpmad"="false" "no-infs-fp-math"="true" "no-nans-fp-math"="true" "no-signed-zeros-fp-math"="true" "no-trapping-math"="true" "probe-stack"="___chkstk_darwin" "stack-protector-buffer-size"="8" "target-cpu"="penryn" "target-features"="+cx16,+cx8,+fxsr,+mmx,+sahf,+sse,+sse2,+sse3,+sse4.1,+ssse3,+x87" "unsafe-fp-math"="true" "use-soft-float"="false" }
attributes #6 = { "correctly-rounded-divide-sqrt-fp-math"="false" "darwin-stkchk-strong-link" "disable-tail-calls"="false" "frame-pointer"="all" "less-precise-fpmad"="false" "no-infs-fp-math"="true" "no-nans-fp-math"="true" "no-signed-zeros-fp-math"="true" "no-trapping-math"="true" "probe-stack"="___chkstk_darwin" "stack-protector-buffer-size"="8" "target-cpu"="penryn" "target-features"="+cx16,+cx8,+fxsr,+mmx,+sahf,+sse,+sse2,+sse3,+sse4.1,+ssse3,+x87" "unsafe-fp-math"="true" "use-soft-float"="false" }
attributes #7 = { argmemonly nofree }
attributes #8 = { nounwind }

!llvm.module.flags = !{!0, !1, !2}
!llvm.ident = !{!3}

!0 = !{i32 2, !"SDK Version", [3 x i32] [i32 10, i32 15, i32 6]}
!1 = !{i32 1, !"wchar_size", i32 4}
!2 = !{i32 7, !"PIC Level", i32 2}
!3 = !{!"Apple clang version 12.0.0 (clang-1200.0.32.29)"}
!4 = !{!5, !5, i64 0}
!5 = !{!"double", !6, i64 0}
!6 = !{!"omnipotent char", !7, i64 0}
!7 = !{!"Simple C++ TBAA"}
!8 = !{!9, !9, i64 0}
!9 = !{!"float", !6, i64 0}
!10 = !{!11}
!11 = distinct !{!11, !12}
!12 = distinct !{!12, !"LVerDomain"}
!13 = !{!14}
!14 = distinct !{!14, !12}
!15 = distinct !{!15, !16}
!16 = !{!"llvm.loop.isvectorized", i32 1}
!17 = distinct !{!17, !16}
!18 = !{!19}
!19 = distinct !{!19, !20}
!20 = distinct !{!20, !"LVerDomain"}
!21 = !{!22, !23, !24}
!22 = distinct !{!22, !20}
!23 = distinct !{!23, !20}
!24 = distinct !{!24, !20}
!25 = !{!22}
!26 = !{!23, !24}
!27 = !{!23}
!28 = !{!24}
!29 = distinct !{!29, !16}
!30 = distinct !{!30, !16}
!31 = !{!32}
!32 = distinct !{!32, !33}
!33 = distinct !{!33, !"LVerDomain"}
!34 = !{!35, !36, !37}
!35 = distinct !{!35, !33}
!36 = distinct !{!36, !33}
!37 = distinct !{!37, !33}
!38 = !{!35}
!39 = !{!36, !37}
!40 = !{!36}
!41 = !{!37}
!42 = distinct !{!42, !16}
!43 = distinct !{!43, !16}
!44 = !{!45}
!45 = distinct !{!45, !46}
!46 = distinct !{!46, !"LVerDomain"}
!47 = !{!48}
!48 = distinct !{!48, !46}
!49 = distinct !{!49, !16}
!50 = distinct !{!50, !16}
!51 = !{!52}
!52 = distinct !{!52, !53}
!53 = distinct !{!53, !"LVerDomain"}
!54 = !{!55}
!55 = distinct !{!55, !53}
!56 = distinct !{!56, !16}
!57 = distinct !{!57, !16}
!58 = !{!59}
!59 = distinct !{!59, !60}
!60 = distinct !{!60, !"LVerDomain"}
!61 = !{!62}
!62 = distinct !{!62, !60}
!63 = distinct !{!63, !16}
!64 = distinct !{!64, !16}
!65 = !{!66}
!66 = distinct !{!66, !67}
!67 = distinct !{!67, !"LVerDomain"}
!68 = !{!69}
!69 = distinct !{!69, !67}
!70 = distinct !{!70, !16}
!71 = distinct !{!71, !16}
!72 = !{!73, !73, i64 0}
!73 = !{!"int", !6, i64 0}
