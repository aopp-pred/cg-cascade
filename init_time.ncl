;===================================================================================================
; uses trusted file to set up 'time_out'
; or vector of desired times
; needed for later when the cloud variables are read in
; and some are missing
;===============================================================

function init_time(generic_path:string,timestep_path:integer,timestep_idx:integer,filename:string,\
                            var_no:string,variable_name:string,flag_interp:logical,t_step_new_ec:integer)

local generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,t_step_new_ec,\
      no_time_slices,time_in,t_count,data_path_matr,data_path,in_file,time_data_tmp,no_time_in,\
      no_steps_tot,no_steps_per_old_step,time_out,i_time,time_tmp,i_st,i_en

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   READ IN GENERIC DATA TO DETERMINE TIME
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; 1. how many time slices do we have in total
  no_time_slices = dimsizes(timestep_idx)
  time_in = new(no_time_slices,integer)

  do t_count = 0,dimsizes(timestep_path)-1

; JASMIN
     data_path_matr = (/generic_path,"/p",tostring(timestep_path(t_count)),"/v0/",var_no,"/",    \
                        filename,"p",tostring(timestep_path(t_count)),"_",var_no,".nc"/)
     data_path      = str_concat(data_path_matr)

  ;  Read in time from file
     in_file        = addfile(data_path,"r")
     time_data_tmp  = in_file->time
     no_time_in = dimsizes(time_data_tmp)

  ; convert CASCADE time from days to seconds
     time_in(t_count) = toint( 60*60*round(time_data_tmp(timestep_idx(t_count))*24,0) )

     if (no_time_in.eq.1) then
        ; time is not defined co-ord variable
        if (timestep_idx(t_count).ne.0) then
           print("##############################################################")
           print("#  Error in setting up time vector                           #")
           print("#  no time in must always be gt 1 for 'trusted' files        #")
           print("##############################################################")
        end if
     end if
  
  end do

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   TEMPORAL INTERPOLATION IF REQUIRED TO SET UP VECTOR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (flag_interp) then
     ; perform linear interpolation in time
     print("    perform temporal interpolation")

     no_steps_tot = round(1+ (time_in(dimsizes(time_in)-1)-time_in(0))  /t_step_new_ec,3)
     no_steps_per_old_step = round(1+(time_in(1)-time_in(0))/t_step_new_ec,3)
     time_out = new(no_steps_tot,typeof(time_in))

     ; ensure time_in times are exactly in time_out
     do i_time=0,dimsizes(time_in)-2   ; minus 1 as start at 0, minus 1 as only interpolate within times
        time_tmp = toint( fspan(time_in(i_time),time_in(i_time+1),no_steps_per_old_step) )
        i_st =    i_time *(no_steps_per_old_step-1)
        i_en = (1+i_time)*(no_steps_per_old_step-1)
        time_out(i_st:i_en) = time_tmp
        delete([/time_tmp,i_st,i_en/])
     end do

  else
     time_out = time_in
  end if

  time_out!0 = "time"
  time_out@calendar = "gregorian"
  time_out@standard_name = "time" ;
  time_out@units = "seconds since 2009-04-06 00:00:00" ;

  return(time_out)

end




