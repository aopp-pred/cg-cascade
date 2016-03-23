; cascade_qtend.ncl

procedure cascade_qtend(fout:file)

local fout, q_mr_data_out, u_data_out, v_data_out, qtend_data_out


begin

  ;==================================================================================
  ; derive advective Q tendency

    print("procedure advective Q tendency")

    ;=====================
    ; load required variables in from file
    q_mr_data_out        = fout->q
    u_data_out           = fout->u
    v_data_out           = fout->v

    ;=====================
    ; calculation

    qtend_data_out = calc_advection(q_mr_data_out,u_data_out,v_data_out)

    ;-- save to file
    add_to_file(fout,qtend_data_out  , "qadv")

    delete([/q_mr_data_out, u_data_out, v_data_out, qtend_data_out/])

end
