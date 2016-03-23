; cascade_ttend.ncl

procedure cascade_ttend(fout:file)

local fout, t_data_out, u_data_out, v_data_out, ttend_data_out


begin

  ;==================================================================================
  ; derive advective T tendency

    print("procedure advective T tendency")

    ;=====================
    ; load required variables in from file
    t_data_out           = fout->t
    u_data_out           = fout->u
    v_data_out           = fout->v

    ;=====================
    ; calculation

    ttend_data_out = calc_advection(t_data_out,u_data_out,v_data_out)

    ;-- save to file
    add_to_file(fout,ttend_data_out  , "tadv")

    delete([/t_data_out, u_data_out, v_data_out, ttend_data_out/])

end
