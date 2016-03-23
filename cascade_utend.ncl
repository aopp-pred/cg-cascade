; cascade_utend.ncl

procedure cascade_utend(fout:file)

local fout, u_data_out, v_data_out, utend_data_out


begin

  ;==================================================================================
  ; derive advective U tendency

    print("procedure advective U tendency")

    ;=====================
    ; load required variables in from file
    u_data_out           = fout->u
    v_data_out           = fout->v

    ;=====================
    ; calculation

    utend_data_out = calc_advection(u_data_out,u_data_out,v_data_out)

    ;-- save to file
    add_to_file(fout,utend_data_out  , "uadv")

    delete([/u_data_out, v_data_out, utend_data_out/])

end
