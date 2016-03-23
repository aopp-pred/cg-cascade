; cascade_vtend.ncl

procedure cascade_vtend(fout:file)

local fout, u_data_out, v_data_out, vtend_data_out


begin

  ;==================================================================================
  ; derive advective V tendency

    print("procedure advective V tendency")

    ;=====================
    ; load required variables in from file
    u_data_out           = fout->u
    v_data_out           = fout->v

    ;=====================
    ; calculation

    vtend_data_out = calc_advection(v_data_out,u_data_out,v_data_out)

    ;-- save to file
    add_to_file(fout,vtend_data_out  , "vadv")

    delete([/u_data_out, v_data_out, vtend_data_out/])

end
