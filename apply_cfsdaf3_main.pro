


pro  apply_cfsdaf3_main
 
  compile_opt idl2
  envi,/restore_base_save_files
  envi_batch_init
  t0=systime(1)

  ;设置预测时刻的存储地
  fine_pre_cfsdaf = "D:\OneDrive\Phd\paper\firstpaper\MFSDAFfusion\testdata\adjustclipbeijing\Default\cfsdaf_12_c_idw_new"                               ;设置最终的预测结果路径
  
  fine_base_t1 = "D:\OneDrive\Phd\paper\firstpaper\MFSDAFfusion\testdata\adjustclipbeijing\Default\adjust_landsat.lst.20000914.dat"                               ;输入t1时刻的fine数据
  coarse_base_t1 = "D:\OneDrive\Phd\paper\firstpaper\MFSDAFfusion\testdata\adjustclipbeijing\Default\30.modis.20000914.dat"                                 ;输入t1时刻的coarse数据
  fine_visdata_t1 = "D:\OneDrive\Phd\paper\firstpaper\MFSDAFfusion\testdata\adjustclipbeijing\Default\adjust_landsat.visdata.20000914.dat"                                ;输入t1时刻的可见光数据
  
  ;端元波普数据
  endmember = "D:\OneDrive\Phd\paper\firstpaper\MFSDAFfusion\testdata\adjustclipbeijing\Default\SVD_Endmembers.csv"                                ;输入端元波谱数据
  coarse_pre_tp = "D:\OneDrive\Phd\paper\firstpaper\MFSDAFfusion\testdata\adjustclipbeijing\Default\30.modis.20021107.dat"                            ;输入tp时刻的coarse数据
  
  ;可变参数
  ratio =32                                   ;设置modis和landsat像素之间的差异值
  w = 2                                       ;设置解混半径大小，若设置2，则窗口大小为2*2+1 = 5
  w_weight = 25                               ;设置距离加权相似像元的窗口
  num_similar_pixel = 30                      ;设置相似像元大小
  temp_file = "d:\temp"                       ;定义缓存路径，中间路径文件。运行前必须在磁盘目录下设置
  interpolation =1                            ;0代表tps插值，1代表IDW插值(默认)
  
  cfsdaf3_main,ratio,w,w_weight,num_similar_pixel,temp_file,interpolation,fine_base = fine_base_t1,coarse_base = coarse_base_t1,coarse_pre = coarse_pre_tp,fine_visdata = fine_visdata_t1,$
    endmember = endmember , fine_pre_cfsdaf = fine_pre_cfsdaf

    ;打印时间
  print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'
end