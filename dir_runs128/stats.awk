BEGIN{ tot_nonCKPTtime = 0.0 ; i = 0 ; num_nonCKPTsteps = 0; num_steps = 0; tot_time = 0.0; max = -999999999  ; min = 9999999999 }
{
    i ++ 
    if ( i > 0 ) {
       tot_time += $1
       num_steps ++
    }
# skip checkpoint steps which for now is set to 45 and 90 for a 90 step run with
# checkpoint every 45 steps.
    if ( i != 45 && i != 90) {
       num_nonCKPTsteps ++
       tot_nonCKPTtime += $1
       if ( $1 > max ) max = $1
       if ( $1 < min ) min = $1
    }
#   printf("%d, %d, %d, %f, %f, %f\n", i, num_steps, num_nonCKPTsteps, $1, tot_time, tot_nonCKPTtime); 
}
# and last
END{printf("Non-checkpoint steps---\n%10s  %8d\n%10s  %15f\n%10s  %15f\n%10s  %15f\n%10s  %15f\n%10s  %15f\n","steps:",num_nonCKPTsteps,"max:",max,"min:",min,"sum:",tot_nonCKPTtime,"mean:",tot_nonCKPTtime/(num_nonCKPTsteps*1.0),"mean/max:",(tot_nonCKPTtime/(num_nonCKPTsteps*1.0))/max); printf("\nAll time steps ---\n     step: %8d, sum: %15f, mean: %15f\n",num_steps,tot_time,(tot_time/(num_steps*1.0))) }
