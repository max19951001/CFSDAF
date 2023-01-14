


;打开数据函数
pro getdata,imgdata = imgdata,ns = ns,nl = nl,nb = nb,data_type = data_type, filename = filename,map_info = map_info, fid = fid, dims = dims
  envi_open_file,filename,r_fid = fid
  envi_file_query,fid,ns = ns,nl = nl,nb = nb,data_type = data_type
  map_info = envi_get_map_info(fid = fid)
  dims = [-1,0,ns - 1 ,0,nl - 1]
  case data_type of
    1:imgdata = bytarr(ns,nl,nb)    ;  byte  byte
    2:imgdata = intarr(ns,nl,nb)    ;  int  integer
    3:imgdata = lonarr(ns,nl,nb)    ;  long  longword integer
    4:imgdata = fltarr(ns,nl,nb)    ;  float  floating point
    5:imgdata = dblarr(ns,nl,nb)    ;  double  double-precision floating
    6:imgdata = complexarr(ns,nl,nb); complex, single-precision, floating-point
    9:imgdata = dcomplexarr(ns,nl,nb);complex, double-precision, floating-point
    12:imgdata = uintarr(ns,nl,nb)   ; unsigned integer vector or array
    13:imgdata = ulonarr(ns,nl,nb)   ;  unsigned longword integer vector or array
    14:imgdata = lon64arr(ns,nl,nb)   ;a 64-bit integer vector or array
    15:imgdata = ulon64arr(ns,nl,nb)   ;an unsigned 64-bit integer vector or array
  endcase
  for i = 0,nb-1 do begin
    dt = envi_get_data(fid = fid,dims = dims,pos=i)
    imgdata[*,*,i] = dt[*,*]
  endfor
end
;*****************************************************************************************************
;用于时间预测的函数
function cost_fun, x
  common v_pub1, ind_v ;ind_v表示丰度值
  common v_pub2, dep_v ;此处dep_v指的是modis像素的变化量
  l=total((ind_v##x - dep_v)^2)
  return, l
end

;用于结合时间和空间预测的解混
function square_weight, x
  common v_pub3, model
  common v_pub4, obs

  g = dblarr(2)
  g[0] = x[0] + x[1]
  g[1] = total((model#x - transpose(obs))^2)
  return, g
end
;****************************************************main program**********************************************************

pro cfsdaf3_main,ratio,w,w_weight,num_similar_pixel,temp_file,interpolation, fine_base = fine_base,coarse_base = coarse_base,coarse_pre = coarse_pre,fine_visdata = fine_visdata,$ 
 endmember= endmember ,  fine_pre_cfsdaf = fine_pre_cfsdaf
 
 
  compile_opt IDL2
  envi, /restore_base_save_files
  envi_batch_init

  ;关键字参数分别为fine_base,coarse_base,coarse_pre,fine_visdata, em, fine_pre,fine_pre_cfsdaf
  ;fine_base,fine_pre,coarse_base,coarse_pre, fine_visdata,em,fine_pre_cfsdaf
  ;
  ;获取到运行文件的目录
  routine_dir=file_dirname(routine_filepath('cfsdaf3_main'))+'\'
  print,routine_dir
  

  ;调用函数设置
  ;***********************************************************************
  ;(1)把数据读进数组中,分别是fine_lst_base,coarse_lst_base,coarse_lst_pre,visdata
  ;读取filename_fine_base,存到fine_lst_base数组中
  getdata,imgdata = fine_lst_base, filename = fine_base,ns = ns_f, nl = nl_f,nb = nb_f, map_info = map_info_f
  getdata,imgdata = coarse_lst_base, filename = coarse_base
  getdata,imgdata = coarse_lst_pre, filename = coarse_pre
  getdata,imgdata = visdata, filename = fine_visdata,fid = fid_vis, ns = ns_vis, nl = nl_vis, nb = nb_vis, map_info = map_info_vis
  ;可见光数据的信息
  ;***************************************************************************
  ;(2)调整异常值,用两倍标准差进行处理
  fine_lst_base_filter = fltarr(ns_f,nl_f,nb_f)
  coarse_lst_base_filter = fltarr(ns_f,nl_f,nb_f)
  coarse_lst_pre_filter = fltarr(ns_f,nl_f,nb_f)
  fine_lst_base_dlt = fltarr(ns_f,nl_f,nb_f)
  coarse_lst_base_dlt = fltarr(ns_f,nl_f,nb_f)
  coarse_lst_pre_dlt = fltarr(ns_f,nl_f,nb_f)
  for ib = 0 , nb_f -1 do begin
    
    ;计算标准差
    fine_lst_base_std= stddev(fine_lst_base[*,*,ib])
    coarse_lst_base_std= stddev(coarse_lst_base[*,*,ib])
    coarse_lst_pre_std = stddev(coarse_lst_pre[*,*,ib])
    
    ;计算平均值,5是邻域半径
    fine_lst_base_filter[*,*,ib] = mean_filter(fine_lst_base[*,*,ib], 5) 
    coarse_lst_base_filter[*,*,ib] = mean_filter(coarse_lst_base[*,*,ib], 5)
    coarse_lst_pre_filter[*,*,ib] = mean_filter(coarse_lst_pre[*,*,ib], 5)
    
    ;计算差值
    fine_lst_base_dlt[*,*,ib] = abs(fine_lst_base[*,*,ib] - fine_lst_base_filter[*,*,ib])
    coarse_lst_base_dlt[*,*,ib] = abs(coarse_lst_base[*,*,ib]-coarse_lst_base_filter[*,*,ib])
    coarse_lst_pre_dlt[*,*,ib] = abs(coarse_lst_pre[*,*,ib]-coarse_lst_pre_filter[*,*,ib])
    
    ww = 3; 设置窗口半径
    for i=0, nl_f -1 do begin
      for j=0, ns_f -1 do begin
        ;
       
        ai = max([0,i-ww])
        bi = min([nl_f-1,i+ww])
        aj = max([0,j-ww])
        bj = min([ns_f-1,j+ww])
        win_size_win = (bj-aj+1)*(bi-ai+1);窗口大小
        ;针对fine_lst_base操作
        if fine_lst_base_dlt[j,i,ib] gt 2*fine_lst_base_std then begin 
          ;获取到窗口内不满足条件的数组的下表，index是不满足条件的数目
          index_yes = where(fine_lst_base_dlt[ai:bi,aj:bj,ib] gt 2*fine_lst_base_std,complement = index, ncomplement = index_value)
          ;将不符合的数据剔除,对剩下的数据进行平均计算
          
          fine_lst_base[j,i,ib] = total((fine_lst_base[ai:bi,aj:bj,ib])[index])/index_value     
        endif
        
        ;针对coarse_lst_base操作
        if coarse_lst_base_dlt[j,i,ib] gt 2*coarse_lst_base_std then begin
        ;获取到窗口内不满足条件的数组的下表
        index_yes = where(coarse_lst_base_dlt[ai:bi,aj:bj,ib] gt 2*coarse_lst_base_std,complement = index, ncomplement = index_value)
        ;将不符合的数据剔除,平均下来计算
        coarse_lst_base[j,i,ib] = total((coarse_lst_base[ai:bi,aj:bj,ib])[index])/index_value
        endif
        
        ;针对coarse_lst_pre操作
        if coarse_lst_pre_dlt[j,i,ib] gt 2*coarse_lst_pre_std then begin
          ;获取到窗口内不满足条件的数组的下表
          index_yes = where(coarse_lst_pre_dlt[ai:bi,aj:bj,ib] gt 2*coarse_lst_pre_std,complement = index, ncomplement = index_value)
          ;将不符合的数据剔除,平均下来计算
          coarse_lst_pre[j,i,ib] = total((coarse_lst_pre[ai:bi,aj:bj,ib])[index])/index_value
        endif
        
   
      endfor
    endfor
  endfor
  ;
  ;(3).光谱解混以及求得modis丰度值
  print,endmember
  em_spec=read_csv(endmember,header=em_name);读取端元
  em_samples=n_elements(em_name) ;返回列数，也就是端元数
  em_lines  =n_elements(em_spec.(0));em_lines等于波段数
  temp=fltarr(em_samples,em_lines)
  for i=0,em_samples-1 do temp[i,*]=float(em_spec.(i))
  em_spec=temporary(temp);em_spec存储各个波段光谱值

  ;fcls,约束最小二乘解混，调用外部程序
  cd,routine_dir
  fbaseabd_fp=routine_dir+'vis'+'_abundance.tif' ;保存解混的风度图
  print,fbaseabd_fp
  ;判断路径中是否存在文件fbaseabd_fp
  result = file_test(fbaseabd_fp)
  if result eq 1 then begin
     fbaseabd_fp=routine_dir+'vis'+'_abundance1.tif'
  endif
  ;cmdstr='abundancecaculatemodule.exe '+adjust_visdata+' '+endmember+' '+fbaseabd_fp ;格式化
  cmdstr=["abundancecaculatemodule.exe",fine_visdata,endmember,fbaseabd_fp];adjust_visdata为临时存储的改正异常值后的值
  print,cmdstr
  spawn,cmdstr,/hide ;spawn调用外部程序，会出现一个窗口进行执行
  envi_open_file,fbaseabd_fp,r_fid=fid_fabd ;fid =fabd_fid,接下来执行就是对fabd_fid执行
  if fid_fabd eq -1 then begin
    envi_batch_exit
    print,'光谱解混失败'
    return
  endif
  print,'最小二乘解混结束' ;屏幕输出这个表示成功解混，得到每个landsat像素内端元的丰度值


  ;求取modis像素的丰度，也就是分辨率960米下的丰度，ratio = 32
  ;首先对fabd_fid执行操作，获取它的参数
  envi_file_query,fid_fabd,ns = ns_fabd, nl = nl_fabd, nb = nb_fabd, dims = dims_fabd
  fabd_img = fltarr(ns_fabd,nl_fabd,nb_fabd)
  ;将landsat丰度值读取到数组fabd_img
  for i = 0,nb_fabd-1 do begin
    dt = envi_get_data(fid = fid_fabd,dims = dims_fabd,pos=i)
    fabd_img[*,*,i] = dt[*,*]
  endfor

  ;ns_c和nl_c分别是粗分辨率下的行列数
  ns_c = ceil(float(ns_f)/ratio);此处向上取整，保证所有的像元都可以遍历到
  nl_c = ceil(float(nl_f)/ratio);获取modis分辨率下的行列号，目的是求这个下的丰度值
  print,nl_c

  ;定义一个数组来存储modis像素丰度
  cabd_img = fltarr(ns_c,nl_c,nb_fabd)
  ;采用循环求取modis丰度值
  for ci = 0, nl_c-1 do begin
    for cj = 0, ns_c-1 do begin
      win_data=fabd_img[(cj*ratio):((cj+1)*ratio-1),(ci*ratio):((ci+1)*ratio-1),*] ;win_data存储的一个MODIS像素内所有Landsat像素的丰度值
      cabd=fltarr(em_samples);cabd存储的临时粗像素的丰度，也就是从landsat聚集到modis中
      for ib=0,em_samples-1 do begin
        cabd[ib]=mean(win_data[*,*,ib])
      endfor
      if total(cabd) lt 0.99 then begin ;将丰度小于1的残差赋值给端元丰度最大的那个
        redi=1.0-total(cabd);残差值
        sortCabd=sort(cabd)
        maxAbundance=cabd[sortCabd[fabd_nb-1]]
        maxAbundance+=redi
        cabd[sortCabd[fabd_nb-1]]=maxAbundance
      endif
      cabd_img[cj,ci,*]=cabd
    endfor
  endfor
  print,'丰度聚集结束'

  ;******************************************************************************************************
  ;(4)进行时间预测的约束解混

  ;定义coarse和fine分辨率下的索引，方便在解混的时候用
  ii=0
  index_f=intarr(ns_f,nl_f)
  index_c=intarr(ns_c,nl_c)
  for i=0, ns_c-1, 1 do begin
    for j=0,nl_c-1,1 do begin
      index_f[i*ratio:(i+1)*ratio-1, j*ratio:(j+1)*ratio-1]=ii
      index_c[i,j]=ii
      ii=ii+1.0
    endfor
  endfor

  ;landsat行列索引值
  row_ind=intarr(ns_f,nl_f)
  col_ind=intarr(ns_f,nl_f)
  for i=0,ns_f-1 do begin
    col_ind[i,*]=i
  endfor
  for i=0,nl_f-1 do begin
    row_ind[*,i]=i
  endfor

  ;获取到聚集后一个modis内变化的土地利用类型占比多少,以及聚集后的粗分辨率数据。定义三个数组
  C_coarse_lst_base = fltarr(ns_c,nl_c,nb_f);存储的基准日期聚集到modis像素中的值
  c_coarse_lst_pre = fltarr(ns_c,nl_c,nb_f);存储的基准日期聚集到modis像素中的值
  c_fine_lst_base = fltarr(ns_c,nl_c,nb_f);对基准日期的数据进行采样
  change_c = fltarr(ns_c,nl_c,nb_f);获取到调整后modis的差值
  row_c=fltarr(ns_c,nl_c)
  col_c=fltarr(ns_c,nl_c)
  for nb =0, nb_f-1 do begin ;考虑多波段
    for ci = 0, nl_c-1 do begin
      for cj = 0, ns_c-1 do begin
        win_data1 = mean(coarse_lst_base[cj*ratio:((cj+1)*ratio-1),ci*ratio:((ci+1)*ratio-1),nb])
        win_data2 = mean(coarse_lst_pre[cj*ratio:((cj+1)*ratio-1),ci*ratio:((ci+1)*ratio-1),nb])
        win_data3 = mean(fine_lst_base[cj*ratio:((cj+1)*ratio-1),ci*ratio:((ci+1)*ratio-1),nb])
        c_coarse_lst_base[cj,ci,nb] = win_data1
        c_coarse_lst_pre[cj,ci,nb] = win_data2
        c_fine_lst_base[cj,ci,nb] = win_data3
;       change_c[cj,ci,nb] = win_data2-win_data1

        ind_c=where(index_f eq index_c[cj,ci])
        row_c[cj,ci]= mean(row_ind[ind_c])
        col_c[cj,ci]= mean(col_ind[ind_c])
      endfor
    endfor
  endfor
  print,"okkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk"

 ;调整传感器引起的差异
 ;存储的是调整后的预测日期的modis数据
 adjust_coarse_lst_pre = fltarr(ns_c,nl_c,nb_f)
 adjust_coarse_lst_base = fltarr(ns_c,nl_c,nb_f)
 
 for ib = 0,nb_f -1 do begin
   y=reform(c_fine_lst_base[*,*,ib],ns_c*nl_c)
   x=reform(c_coarse_lst_base[*,*,ib],ns_c*nl_c)
   coefs=linfit(x,y) ;coefs[1] = a, coefs[0] = b
   adjust_coarse_lst_pre[*,*,ib]= coefs[1]*c_coarse_lst_pre[*,*,ib]+coefs[0]
   adjust_coarse_lst_base[*,*,ib] = coefs[1]*c_coarse_lst_base[*,*,ib]+coefs[0]
   print,coefs
 endfor

    
 for ib = 0 ,nb_f -1 do begin
   change_c[*,*,ib] = adjust_coarse_lst_pre[*,*,ib]-adjust_coarse_lst_base[*,*,ib] 
 endfor
 
  ;开始进行时间预测解混，最小二乘线性解混
  ;定义解混的最小最大值,nb_f为波段，em_samples为端元数
  min_allow=fltarr(em_samples,nb_f)
  max_allow=fltarr(em_samples,nb_f)
  ;定义时间变化值
  temporal_change = fltarr(ns_f,nl_f,nb_f)

  common v_pub1
  common v_pub2
  gbnd    =[-50,50]      ;设置两个日期内modis的变化范围
  nobj    = 0           ;因为所求结果为1个，设置索引为0
  lcomp   = 'cost_fun'

  for ci = 0, nl_c -1 do begin
    for cj = 0, ns_c -1 do begin
      ai = max([0,ci-w])
      bi = min([nl_c-1,ci+w])
      aj = max([0,cj-w])
      bj = min([ns_c-1,cj+w])

      search_win = (bj-aj+1)*(bi-ai+1) ;搜索窗口大小
      ;设置modis对应下landsat窗口的位置
      fai=ci*ratio
      fbi=(ci+1)*ratio-1
      faj=cj*ratio
      fbj=(cj+1)*ratio-1

      ;对于边界像素进行解混的话，这个窗口内像素往往是不够的，因此需要进行调整
      total_win = (w*2+1)*(w*2+1);total_win为完整的窗口
      ;对于搜索窗口不够解混的窗口，进行扩展;
      while search_win lt total_win do begin
        if (ci-w)le 0 then bi++
        if (ci+w)ge nl_c then ai--
        if (cj-w)le 0 then bj++
        if (cj+w)ge ns_c then aj--
        search_win = (bj-aj+1)*(bi-ai+1)
      endwhile

      fabd_img_temp = fabd_img[faj:fbj,fai:fbi,*] ;将窗口内的端元丰度存储到临时的fabd_img_temp，为了便于最终的求解变化量
      ind_v=transpose(reform(cabd_img[aj:bj,ai:bi,*],search_win,em_samples));modis像素每一类的丰度值,然后将其转置
      for ib =0, nb_f-1 do begin
        min_allow[*,ib] = min(change_c[aj:bj,ai:bi,ib]-stddev(change_c[aj:bj,ai:bi,ib])) ;允许的粗像素的空间变化范围
        max_allow[*,ib] = max(change_c[aj:bj,ai:bi,ib]+stddev(change_c[aj:bj,ai:bi,ib]))
        dep_v = reform(change_c[aj:bj,ai:bi,ib],search_win) ;设置变化值
        x  = fltarr(1,em_samples);最终的端元变化结果存在这个里面
        xbnd  = [[min_allow[*,ib]],[max_allow[*,ib]]];设置变量的变化范围，也就是每一端元的变化范围
        constrained_min, x, xbnd, gbnd, nobj, lcomp,inform, nstop = 5 ;此处x返回的是每个modis像素内的端元变化范围
        ds_change=fabd_img_temp*rebin(reform(x,1,1,em_samples),ratio,ratio,em_samples,/sample);最邻近采样为了保证值不变的情况下进行矩阵相乘
        temporal_change[faj:fbj,fai:fbi,ib]=total(ds_change,3);
      endfor
    endfor
  endfor
  print,"时间预测结束"
  
;  t0=systime(1)
 
  ;(5)使用tps插值/IDW进行空间预测
  ;定义空间变化量
  spatial_change = fltarr(ns_f,nl_f,nb_f)
  ;根据变量interpolation来控制最终的插值结果，如果interpolation = 0（默认值）,采用tps插值；如果interpolation = 1,则采用idw插值
  if interpolation eq 0 then begin
    for ib=0,nb_f-1 do begin
      print,"ok"
      spatial_change[*,*,ib] = min_curve_surf(change_c[*,*,ib],col_c, row_c,/tps, xpout=col_ind,  ypout=row_ind)
    endfor    
  endif
  
  ;使用IDW插值
  if interpolation eq 1 then begin
    ;定义偏移值，确定每个粗像素在一个细像素中的位置
    bias = floor(ratio/2);
    ;循环
    for ci=0, nl_c -1 do begin
      for cj=0,ns_c -1 do begin
        ;确定每次细像素的起始位置
        fi = w*ci
        fj = w*cj
        ;确定粗像素的窗口大小
        cai = max([0,ci-1])
        cbi = min([nl_c-1,ci+1])
        caj = max([0,cj-1])
        cbj = min([ns_c-1,cj+1])

        ;获取到粗分辨率对应的细分辨率数据的索引值
        fai=ci*ratio
        fbi=(ci+1)*ratio-1
        faj=cj*ratio
        fbj=(cj+1)*ratio-1

        for ib = 0, nb_f-1 do begin

          ;定义窗口内的所有modis像素
          coarse_value = change_c[caj:cbj,cai:cbi,ib] ;将值赋值过来
          coarse_value = reform(coarse_value,(cbi-cai+1)*(cbj-caj+1))

          ;确定窗口内每一粗像素在整个影像细分辨率像素处的索引
          coarse_index = intarr(cbj-caj+1,cbi-cai+1,2)

          for i = 0, cbi-cai do begin
            for j = 0, cbj-caj do begin
              coarse_index[j,i,0] = (ci)*ratio+bias
              coarse_index[j,i,1] = (cj)*ratio+bias
            endfor
          endfor

          coarse_index_adjust = reform(coarse_index,(cbi-cai+1)*(cbj-caj+1),2)
          ;定义存储距离的变量
          distance = fltarr((cbi-cai+1)*(cbj-caj+1))
          ;定义存储反距离平方的变量
          inver_distance = fltarr((cbi-cai+1)*(cbj-caj+1))
          ;定义每一个modis对于最终预测权重值
          weight = fltarr((cbi-cai+1)*(cbj-caj+1))

          ;针对此ci,cj里面的每一个细数据进行求取距离
          for fi = fai, fbi do begin ;fai,fbi分别是一个modis里面待插值的landsat像素位置
            for fj = faj, fbj do begin

              for ib = 0 ,nb_f-1 do begin
                ;temp 为距离差,如果fj,fi刚好在中心modis位置上
                temp = fltarr((cbi-cai+1)*(cbj-caj+1),2)
                temp[*,0] = coarse_index_adjust[*,0]-fi
                temp[*,1] = coarse_index_adjust[*,1]-fj

                for num = 0, (cbi-cai+1)*(cbj-caj+1) -1 do begin
                  distance[num] = sqrt(temp[num,0]*temp[num,0]+temp[num,1]*temp[num,1])  ;按循环求取每一个粗像素到(fj,fi)的距离
                  inver_distance[num] = 1/(distance[num]*distance[num])   ;距离平方的倒数
                endfor

                ;计算反距离的权重
                for num = 0, (cbi-cai+1)*(cbj-caj+1) -1 do begin
                  weight[num] = inver_distance[num]/total(inver_distance)       ;反距离加权
                endfor

                ;计算最终结果，coarse_value*weight
                spatial_change[fj,fi] =transpose(weight)#coarse_value
                ; print,spatial_change[fj,fi]
              endfor
            endfor
          endfor
          ;      print,"一个粗分辨率像素预测结束"
        endfor
      endfor
    endfor

    print,"IDW插值结束"
    ; 但是这一过程导致每个粗像元的中心像素出现了异常值，我们用邻域内像素进行填充
    for ib = 0, nb_f-1 do begin
      result_temp = spatial_change[*,*,ib]
      index = where(finite(result_temp,/nan),count)
      result_temp[index] = 0.0

      for fi = 0, nl_f -1  do begin
        for fj = 0 ,ns_f -1  do begin
          ai = max([0,fi-w])
          bi = min([nl_f-1,fi+w])
          aj = max([0,fj-w])
          bj = min([ns_f-1,fj+w])

          num = (bj-aj+1)*(bi-ai+1)-1 ;窗口内所有的像素，不包括中心像素
          if result_temp[fj,fi] eq 0.0 then begin
            result_temp[fj,fi] = total(result_temp[aj:bj,ai:bi])/num
          endif

        endfor
      endfor
      spatial_change[*,*,ib] = result_temp
    endfor
 
  endif
  print,"空间预测结束"
  
  ;插值结束时间
;  print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'

  spatial_output = routine_dir+ "spatial_change"
  envi_write_envi_file,spatial_change, out_name = spatial_output, ns = ns_f, nl  = nl_f, nb = nb_f, map_info = map_info_f
  

  ;(6)cls解混。
  ;结合时间预测和空间预测一起进行结果
  ;此处的变量是model[0]和model[1]
  
  common v_pub3, model
  common v_pub4, obs ;obs设置目标函数的变化
  fixed_change = fltarr(ns_f,nl_f,nb_f);定义修改后的变化值,在landsat尺度
  
  ;将时间预测和空间预测的结果进行
  temp_spatial_change = fltarr(ns_c,nl_c,nb_f) ;空间预测的聚集
  temp_temporal_change = fltarr(ns_c,nl_c,nb_f) ;时间预测的聚集
  for ci = 0, nl_c-1 do begin
    for cj = 0, ns_c-1 do begin
      ;设置modis对应下landsat窗口的位置
      fai=ci*ratio
      fbi=(ci+1)*ratio-1
      faj=cj*ratio
      fbj=(cj+1)*ratio-1
      
      for ib =0,nb_f -1 do begin       
        temp_spatial_change[cj,ci,ib] = mean(spatial_change[faj:fbj,fai:fbi,ib])
        temp_temporal_change[cj,ci,ib] = mean(temporal_change[faj:fbj,fai:fbi,ib])     
      endfor   
    endfor
  endfor
  
; 
  ;1.结合时间预测和空间预测在modis像素尺度进行
  for ci = 0, nl_c-1 do begin
    for cj = 0, ns_c-1 do begin
      
      ai = max([0,ci-w])
      bi = min([nl_c-1,ci+w])
      aj = max([0,cj-w])
      bj = min([ns_c-1,cj+w])

      search_win = (bj-aj+1)*(bi-ai+1) ;搜索窗口大小
      ;设置modis对应下landsat窗口的位置
      fai=ci*ratio
      fbi=(ci+1)*ratio-1
      faj=cj*ratio
      fbj=(cj+1)*ratio-1

      ;对于边界像素进行解混的话，这个窗口内像素往往是不够的，因此需要进行调整
      total_win = (w*2+1)*(w*2+1);total_win为完整的窗口
      ;对于搜索窗口不够解混的窗口，进行扩展;
      while search_win lt total_win do begin
        if (ci-w)le 0 then bi++
        if (ci+w)ge nl_c then ai--
        if (cj-w)le 0 then bj++
        if (cj+w)ge ns_c then aj--
        search_win = (bj-aj+1)*(bi-ai+1)
      endwhile
      
      
      model = fltarr(search_win,2)
      for ib = 0, nb_f -1 do begin
       
        obs = change_c[aj:bj,ai:bi,ib];在一个modis像素内进行运行
        model[*,0]=reform(temp_temporal_change[aj:bj,ai:bi,ib],search_win)
        model[*,1]=reform(temp_spatial_change[aj:bj,ai:bi,ib],search_win)

        ;定义最小二乘约束解混参数
        xlb=fltarr(2, 1)
        xub=fltarr(2, 1) + 1
        xbnd = [[xlb], [xub]]
        x = [0.5, 0.5] ;初始值
        nobj = 1
        gcomp = 'square_weight'
        gbnd = [[1,0], [1,0]] ;边界忽略，约束g[0]和g[1]
        constrained_min, x, xbnd, gbnd, nobj, gcomp, inform, epstop = 1.0e-8, nstop = 10
        result = x ;将x存到result中，这是一个modis像素下的权重值
        model[*, 0] = model[*, 0] + 0.15*result[0]*stddev(model[*, 0])*randomn(seed, n_elements(obs), 1)
        model[*, 1] = model[*, 1] + 0.15*result[1]*stddev(model[*, 0])*randomn(seed, n_elements(obs), 1)
        constrained_min, x, xbnd, gbnd, nobj, gcomp, inform, epstop = 1.0e-8, nstop = 10
        
        result = x
;        print,x
        ;将一个modis下的变化值赋值给对应的landsat数据中，此刻数据没有考虑传感器的差异
        fixed_change[faj:fbj,fai:fbi,ib]= result[0]*temporal_change[faj:fbj,fai:fbi,ib]+result[1]*spatial_change[faj:fbj,fai:fbi,ib]
      endfor
      
    endfor
  endfor
  
;  ;2,结合时间预测和空间预测在fine像素尺度上进行，这种方法忽略了空间自相关性
;  change_f = coarse_lst_pre - coarse_lst_base ;定义的是coarse像素在fine分辨率的变化值
;  
;  
;  ;在fine像素下，进行每个modis像素拟合
;  for ci = 0, nl_c -1 do begin
;    for cj = 0, ns_c -1 do begin
;      ;设置ci，cj对应下landsat的像素索引
;      fai=ci*ratio
;      fbi=(ci+1)*ratio-1
;      faj=cj*ratio
;      fbj=(cj+1)*ratio-1
;
;      ;定义解混的时间和空间变化量
;      model = fltarr(ratio*ratio,2)
;      for ib = 0, nb_f -1 do begin
;        obs = change_f[faj:fbj,fai:fbi,ib];在一个modis像素内进行运行
;        model[*,0]=reform(temporal_change[faj:fbj,fai:fbi,ib],ratio*ratio)
;        model[*,1]=reform(spatial_change[faj:fbj,fai:fbi,ib],ratio*ratio)
;
;        ;定义最小二乘约束解混参数
;        xlb=fltarr(2, 1)
;        xub=fltarr(2, 1) + 1
;        xbnd = [[xlb], [xub]]
;        x = [0.5, 0.5] ;初始值
;        nobj = 1
;        gcomp = 'square_weight'
;        gbnd = [[1,0], [1,0]] ;边界忽略，约束g[0]和g[1]
;        constrained_min, x, xbnd, gbnd, nobj, gcomp, inform, epstop = 1.0e-8, nstop = 10        
;        result = x ;将x存到result中，这是一个modis像素下的权重值
;        model[*, 0] = model[*, 0] + 0.15*result[0]*stddev(model[*, 0])*randomn(seed, n_elements(obs), 1)
;        model[*, 1] = model[*, 1] + 0.15*result[1]*stddev(model[*, 0])*randomn(seed, n_elements(obs), 1)
;        constrained_min, x, xbnd, gbnd, nobj, gcomp, inform, epstop = 1.0e-8, nstop = 10
;        
;        ;将一个modis下的变化值赋值给对应的landsat数据中，此刻数据没有考虑传感器的差异
;        fixed_change[faj:fbj,fai:fbi,ib]= result[0]*temporal_change[faj:fbj,fai:fbi,ib]+result[1]*spatial_change[faj:fbj,fai:fbi,ib]
;      endfor
;    endfor
;  endfor

  ;存储变化值
  fixed_change_result = routine_dir+"fixed_change_result5"
  envi_write_envi_file,fixed_change,out_name = fixed_change_result,ns = ns_f, nl  = nl_f, nb = nb_f, map_info = map_info_f

  print,"时间预测和空间预测结合结束"
  
  ;定义最终的变化
  last_change = fltarr(ns_f,nl_f,nb_f)
  ;残差纠正，考虑到了我们在细分辨率下解混完全是为了保证计算的效率
  for ci = 0, nl_c -1 do begin
    for cj = 0, ns_c -1 do begin
      ;设置ci，cj对应下landsat的像素索引
      fai=ci*ratio
      fbi=(ci+1)*ratio-1
      faj=cj*ratio
      fbj=(cj+1)*ratio-1

      for ib = 0, nb_f -1 do begin
        temp = change_c[cj,ci,ib]- mean(fixed_change[faj:fbj,fai:fbi,ib]) ;temp存储的每个像元的残差
        last_change[faj:fbj,fai:fbi,ib] = fixed_change[faj:fbj,fai:fbi,ib]+rebin(reform(temp,1,1),ratio,ratio)
      endfor
    endfor
  endfor
  
  last_change_result = routine_dir+"last_change_result5"
  envi_write_envi_file,last_change,out_name = last_change_result,ns = ns_f, nl  = nl_f, nb = nb_f, map_info = map_info_f
  ;*****************************************************************************************************************************
 
  ;(7)采用starfm或者fsdaf的距离加权进行预测

  ;最终的权重加权
  ;设置输出结果
  fine_lst_pre = fltarr(ns_f,nl_f,nb_f)
  ;d_d_all参考FSDAF中的计算

  d_d_all=((w_weight-indgen(w_weight*2+1)#(intarr(1,w_weight*2+1)+1))^2+(w_weight-(intarr(w_weight*2+1)+1)#indgen(1,w_weight*2+1))^2)^0.5
  d_d_all=reform(d_d_all,(w_weight*2+1)*(w_weight*2+1))

  ;按照热红外数据设置寻找相似像元的条件
  similar_th=fltarr(nb_f)
  for iband=0,nb_f-1 do begin
    similar_th[iband]=stddev(fine_lst_base[*,*,iband])*2.0/em_samples ;FSDAF可能直接用地表温度数据来进行搜索相似像素了
  endfor
;  ;按照可见光数据搜索相似像元
;  similar_th = fltarr(nb_vis)
;  for iband = 0, nb_vis-1 do begin
;    similar_th[iband] = stddev(visdata[*,*,iband])*2.0/em_samples  ;此处用可见光数据搜索相似像元
;  endfor
  ;开始了
  for fi = 0, nl_f -1 do begin
    for fj = 0, ns_f -1 do begin
      ;设置窗口大小
      ai=max([0,fi-w_weight])
      bi=min([nl_f-1,fi+w_weight])
      aj=max([0,fj-w_weight])
      bj=min([ns_f-1,fj+w_weight])

      col_wind = indgen(bj-aj+1)#(intarr(1,bi-ai+1)+1)*1.0 ;列索引，就是每列都一样的
      row_wind = (intarr(bj-aj+1)+1)#indgen(1,bi-ai+1)*1.0

      ;寻找窗口内相似像素
      similar_cand = fltarr((bj-aj+1)*(bi-ai+1));放置目标像元与窗口内相似程度,默认为0.0
      position_cand = intarr((bj-aj+1)*(bi-ai+1))+1 ;放置每个相似像元的位置，默认从1开始
      ;按照地表温度数据搜索
      for ib = 0, nb_f-1 do begin
        cand_band = intarr((bj-aj+1)*(bi-ai+1))
        wind_fine = fine_lst_base[aj:bj,ai:bi,ib] ;窗口内的基准日期像素
        s_s = abs(wind_fine-fine_lst_base[fj,fi,ib]) ;计算窗口内像素到中心像素的光谱差，此处是地温差
        similar_cand = similar_cand+s_s/(fine_lst_base[fj,fi,ib]+0.00000001);存储的是相似程度，值越大，越不相似,后面需要从小到大排序
        ind_cand = where(s_s lt similar_th[ib]);判断窗口的光谱差小于给定阈值的索引
        cand_band[ind_cand] = 1;将相似的值赋值1
        position_cand = position_cand*cand_band;position_cand中与中心像元相似的像素赋值为1
      endfor     
      ;
;      ;利用可见光数据进行搜索相似 像元
;      for ib = 0, nb_vis-1 do begin
;        cand_band = intarr((bj-aj+1)*(bi-ai+1))
;        wind_fine = visdata[aj:bj,ai:bi,ib] ;窗口内的基准日期像素
;        s_s = abs(wind_fine-visdata[fj,fi,ib]) ;计算窗口内像素到中心像素的光谱差，此处是地温差
;        similar_cand = similar_cand+s_s/(visdata[fj,fi,ib]+0.00000001);存储的是相似程度，值越大，越不相似,后面需要从小到大排序
;        ind_cand = where(s_s lt similar_th[ib]);判断窗口的光谱差小于给定阈值的索引
;        cand_band[ind_cand] = 1;将相似的值赋值1
;        position_cand = position_cand*cand_band;position_cand中与中心像元相似的像素赋值为1
;      endfor

      indcand = where(position_cand ne 0, number_cand0) ;筛选相似像素的索引，number-cand0为相似像素的个数
      order_dis = sort(similar_cand[indcand]);order_dis存储的是相似像素的从小到大索引值
      number_cand = min([number_cand0,num_similar_pixel]) ;进行筛选
      ind_same_class = indcand[order_dis[0:number_cand-1]];筛选出最相似的像元,为窗口内的索引值


      ;计算相似像元到中心像素的距离
      if((bj-aj+1)*(bi-ai+1) lt (w_weight*2.0+1)*(w_weight*2.0+1)) then begin ;对于图像边缘不够一个窗口的处理
        d_d_cand = ((fj-col_wind[ind_same_class])^2+(fi-row_wind[ind_same_class])^2)^0.5+0.000001
      endif else begin
        d_d_cand = d_d_all[ind_same_class]
      endelse

      ;归一化距离
      d_d_cand=(1.0+d_d_cand/w_weight)*(similar_cand[ind_same_class]+1.0)
      c_d=1.0/d_d_cand
      weight=c_d/total(c_d) ;weight为权重

      for ib = 0 ,nb_f -1 do begin
        fixed_change_temp = ( last_change[aj:bj,ai:bi,ib])[ind_same_class]
        fine_lst_pre[fj,fi,ib] = fine_lst_base[fj,fi,ib]+total(weight*fixed_change_temp)
      endfor

    endfor
  endfor
  print,"距离加权运算结束"
  ;*************************************************************************************************************************
  ;对最终的预测结果进行处理
  
  fine_lst_pre_filter = fltarr(ns_f,nl_f,nb_f)
  fine_lst_pre_dlt = fltarr(ns_f,nl_f,nb_f)
  
  for ib = 0 , nb_f -1 do begin

    ;计算标准差
    fine_lst_pre_std= stddev(fine_lst_pre[*,*,ib])
   

    ;计算平均值,5是邻域半径
    fine_lst_pre_filter[*,*,ib] = mean_filter(fine_lst_pre[*,*,ib], 5)
    

    ;计算差值
    fine_lst_pre_dlt[*,*,ib] = abs(fine_lst_pre[*,*,ib] - fine_lst_pre_filter[*,*,ib])
   
    ww = 3; 设置窗口半径
    for i=0, nl_f -1 do begin
      for j=0, ns_f -1 do begin
        ;
        ai = max([0,i-ww])
        bi = min([nl_f-1,i+ww])
        aj = max([0,j-ww])
        bj = min([ns_f-1,j+ww])
        win_size_win = (bj-aj+1)*(bi-ai+1);窗口大小
        ;针对fine_lst_base操作
        if fine_lst_pre_dlt[j,i,ib] gt 2*fine_lst_pre_std then begin
          ;获取到窗口内不满足条件的数组的下表，index是不满足条件的数目
          index_yes = where(fine_lst_pre_dlt[ai:bi,aj:bj,ib] gt 2*fine_lst_pre_std,complement = index, ncomplement = index_value)
          print,"complements= ",index
          ;将不符合的数据剔除,对剩下的数据进行平均计算
          fine_lst_pre[j,i,ib] = total((fine_lst_pre[ai:bi,aj:bj,ib])[index])/index_value
        endif
       endfor
     endfor
   endfor
  ;print,fine_lst_pre
  envi_write_envi_file,fine_lst_pre, out_name = fine_pre_cfsdaf, map_info = map_info_f , ns =ns_f, nl = nl_f, nb = nb_f

end