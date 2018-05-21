; cascade_utend.ncl

;================================================================================
;   Copyright 2018 Hannah M. Christensen
;
;   Licensed under the Apache License, Version 2.0 (the "License");
;   you may not use this file except in compliance with the License.
;   You may obtain a copy of the License at
;
;       http://www.apache.org/licenses/LICENSE-2.0
;
;   Unless required by applicable law or agreed to in writing, software
;   distributed under the License is distributed on an "AS IS" BASIS,
;   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
;   See the License for the specific language governing permissions and
;   limitations under the License.
;================================================================================


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

    utend_data_out@long_name = "Advective U tendency" ;
    utend_data_out@units     = "m/s^2" ;

    add_to_file(fout,utend_data_out  , "uadv")

    delete([/u_data_out, v_data_out, utend_data_out/])

end
