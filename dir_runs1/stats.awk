BEGIN{ tot_nonIOtime = 0.0 ; i = 0 ; num_nonIOsteps = 0; num_steps = 0; tot_time = 0.0; max = -999999999  ; min = 9999999999 }
{
    i ++ 
    if ( i > 0 ) {
       tot_time += $1
       num_steps ++
    }
# skip IO/checkpoint steps, which for now is set to 10 and 20 for a 20 step run with
# IO/checkpoint every 10 steps.
    if ( i != 10 && i != 20) {
       num_nonIOsteps ++
       tot_nonIOtime += $1
       if ( $1 > max ) max = $1
       if ( $1 < min ) min = $1
    }
#   printf("%d, %d, %d, %f, %f, %f\n", i, num_steps, num_nonIOsteps, $1, tot_time, tot_nonIOtime); 
}
# and last
END{printf("Non-I/O steps---\n%10s  %8d\n%10s  %15f\n%10s  %15f\n%10s  %15f\n%10s  %15f\n%10s  %15f\n","steps:",num_nonIOsteps,"max:",max,"min:",min,"sum:",tot_nonIOtime,"mean:",tot_nonIOtime/(num_nonIOsteps*1.0),"mean/max:",(tot_nonIOtime/(num_nonIOsteps*1.0))/max); printf("\nAll time steps ---\n     step: %8d, sum: %15f, mean: %15f\n",num_steps,tot_time,(tot_time/(num_steps*1.0))) }
