;===================================================================================================
; read in atmospheric data from CASCADE files
;
; can read in from multiple time steps and place in same array
; copy coordinate variables from CASCADE array to finished array
;
; performs temporal interpolation if flag_interp = .TRUE.
;===============================================================

function read_cascade(generic_path:string,timestep_path:string,filename:string,\
                            var_no:string,variable_name:string,flag_interp:logical,t_step_new_ec:float)

local generic_path,timestep_path,filename,var_no,variable_name,t_count,data_path_matr,\
      data_path,in_file,in_data_tmp,data_size,no_times_in_file,in_data,t_st,t_en,var_dims,\
      i_dim,time_tmp,time_flag,time_flag_save,t_en_old,no_time_slices,time_in,time_out,\
      time_dim,no_steps_tot,no_steps_per_old_step,i_time,time_tmp,i_st,i_en,time_days


begin
;  setfileoption("nc","Format","NetCDF4Classic")
;  setfileoption("nc","Compression",1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   READ IN DATA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; 1. work out how many time slices we have in total
  no_time_slices = 0
  do t_count = 0,dimsizes(timestep_path)-1

     data_path_matr = (/generic_path,"/",timestep_path(t_count),"/",filename,timestep_path(t_count),"_",var_no,".nc"/)
     data_path      = str_concat(data_path_matr)

  ;  Input variable dimensions: (time, lat, lon) ;

     in_file = addfile(data_path,"r")
     time_data_tmp    = in_file->time
     time_size = dimsizes(time_data_tmp)
     no_time_slices = no_time_slices + time_size

     delete([/data_path_matr,data_path,in_file,time_data_tmp,time_size/])
   end do

  t_en_old = -1
  time_flag_save = "false"
  do t_count = 0,dimsizes(timestep_path)-1

     data_path_matr = (/generic_path,"/",timestep_path(t_count),"/",filename,timestep_path(t_count),"_",var_no,".nc"/)
     data_path      = str_concat(data_path_matr)
  
  ;  Input variable dimensions: (time, latitude, longitude) ;
     in_file = addfile(data_path,"r")
     in_data_tmp    = in_file->$variable_name$
  
     data_size = dimsizes(in_data_tmp)

     if (dimsizes(data_size).eq.3) then           ; time,lat, lon
        in_data_tmp!1 = "lat"
        in_data_tmp!2 = "lon"
     else
        if (dimsizes(data_size).eq.4) then        ; time, lev, lat, lon
           in_data_tmp!2 = "lat"
           in_data_tmp!3 = "lon"
        end if
     end if



     if iscoord(in_data_tmp,"time")
        ; time is defined co-ord variable
        no_times_in_file = dimsizes(in_data_tmp&time)
     else
        ; time is not defined co-ord variable
        no_times_in_file = 1
     end if

     ; determine if this file has one or more time slices in
     if (data_size(0).eq.no_times_in_file) then   ; orig file has more than one time slice in
        time_flag = "false"
     else                                         ; orig file has just one time in
        time_flag      = "true"
     end if

     ; set up the new empty array into which we will put our multiple time data
     if (t_count.eq.0) then
        var_dims = getvardims(in_data_tmp)
        if time_flag.eq."false" then   ; orig file has more than one time slice in
           data_size(0) = no_time_slices       ; OLD ; data_size(0)*dimsizes(timestep_path)
           in_data = new(data_size,double)
           do i_dim = 0,dimsizes(var_dims)-1
              in_data!i_dim = in_data_tmp!i_dim
           end do
           time_tmp = new(no_time_slices,double)
        else                                         ; orig file has just one time in
           data_size_new = array_append_record(no_time_slices,data_size,0)
;           data_size_new = array_append_record(dimsizes(timestep_path),data_size,0)
           delete(data_size)
           data_size = data_size_new
           delete(data_size_new) 
           in_data = new(data_size,double)
           in_data!0 = "time"
           do i_dim = 0,dimsizes(var_dims)-1
              in_data!(i_dim+1) = in_data_tmp!(i_dim)
           end do
         ;;;;  time_tmp = new(dimsizes(timestep_path),double)
           time_tmp = new(no_time_slices,double)
        end if
        delete([/var_dims,i_dim/])
     end if

     ;- if orig array has only one time in, need to introduce time dimension
     if (time_flag.eq."true")
;        printVarSummary(in_data_tmp)

        dsizes = dimsizes(in_data_tmp)
        dsizes2 = array_append_record(dsizes,1,0)
        tmparray = new(dsizes2,typeof(in_data_tmp))
        data_tmp = conform(tmparray,in_data_tmp,(/0,1,2/))
        ; consider different cases
        if (dimsizes(dsizes).eq.2) then ;lat lon only
           data_tmp!0 = in_data_tmp!0
           data_tmp!1 = in_data_tmp!1
           data_tmp!2 = "time"
           data_tmp2 = data_tmp(time | :, lat | :, lon | :)
           delete(in_data_tmp)
           in_data_tmp = data_tmp2
           delete(data_tmp2)
        else
           if (dimsizes(dsizes).eq.3) then ; lev, lat lon
              data_tmp!0 = in_data_tmp!0
              data_tmp!1 = in_data_tmp!1
              data_tmp!2 = in_data_tmp!2
              data_tmp!3 = "time"
              data_tmp2 = data_tmp(time | :, eta_hh | :, lat | :, lon | :)
              delete(in_data_tmp)
              in_data_tmp = data_tmp2
              delete(data_tmp2)
           end if
        end if
        delete([/dsizes,dsizes2,tmparray,data_tmp/])
     end if

     t_st     = t_en_old+1
     t_en     = t_st+no_times_in_file-1
     t_en_old = t_en

;     t_st     =     t_count*no_times_in_file
;     t_en     = (t_count+1)*no_times_in_file - 1
  
     if (dimsizes(data_size).eq.3) then           ; time,      lat, lon
        in_data(t_st:t_en,:,:) = in_data_tmp
        in_data!1 = "lat"
        in_data!2 = "lon"
     else
        if (dimsizes(data_size).eq.4) then        ; time, lev, lat, lon
           in_data(t_st:t_en,:,:,:) = in_data_tmp
           in_data!2 = "lat"
           in_data!3 = "lon"
        end if
     end if

     ; read in times in days - convert to seconds
     time_days            = in_file->time
     time_tmp(t_st:t_en)  = round(time_days*24*60*60,0)

     delete([/data_path_matr,data_path,in_file,in_data_tmp,no_times_in_file,t_st,t_en,data_size,time_flag,time_days/])
  
  end do

  in_data&time = time_tmp
 
  delete([/time_tmp,t_en_old/])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   TEMPORAL INTERPOLATION IF REQUIRED
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (flag_interp) then
     ; perform linear interpolation in time
     print("    perform temporal interpolation")

     time_in = in_data&time  ; in seconds
     no_steps_tot = round(1+ (time_in(dimsizes(time_in)-1)-time_in(0))  /t_step_new_ec,3)
     no_steps_per_old_step = round(1+(time_in(1)-time_in(0))/t_step_new_ec,3)

     time_out = new(no_steps_tot,typeof(time_in))

     ; ensure time_in times are exactly in time_out
     do i_time=0,dimsizes(time_in)-2   ; minus 1 as start at 0, minus 1 as only interpolate within times
        time_tmp = fspan(time_in(i_time),time_in(i_time+1),no_steps_per_old_step)
        i_st =    i_time *(no_steps_per_old_step-1)
        i_en = (1+i_time)*(no_steps_per_old_step-1)
        time_out(i_st:i_en) = time_tmp
        delete([/time_tmp,i_st,i_en/])
     end do

     time_dim = 0
     in_data_ti = linint1_n_Wrap (in_data&time,in_data,False,time_out,0,time_dim)
     in_data_ti!0 = "time"

     delete([/time_in,no_steps_tot,time_out,in_data/])
  else
     in_data_ti = in_data
     delete([/in_data/])
  end if


  return(in_data_ti)

end




