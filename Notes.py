My thoughts:
    1. divide the variance by the sqrt of the number stacked and
       the mean by the number stacked.
           Issue: The variance didn't fall quickly enough. And in any
                  case the stacking of the clusters should have increased
                  according to sqrt n anyways

    2. create just a signal map and just a noise map. divide each by its
       factor 1/n and 1/sqrt(n). 
            Issue: This won't work, essentially we are scaling down both
                   the noise and the signal and we end up with the same as
                   we would have gotten if we had just stacked one cluster

    3. create just a signal map and just a noise map. stack them and see
       what happens. noise should increase at lower rate than signal and
       in general variance should fall
            Issue: the average variance does rise, but the signal doesn't keep
                   up for some reason

    4. create just a signal map and just a noise map. divide the noise by
       1/ sqrt(n) factor. noise should converge to the variance already set
       while the signal just keeps getting bigger
           Issue: seems to be giving a reasonable result, however we are not
                  seeing a consistant rise in the SNR

Converting Mpc (R500) to arcmins

To-do:
    for the SNR maps:
        -test with just noise to attempt to detect possible bug
        -test with different size bins
        -test with different resoultion maps, smaller or larger pixel size
    for different types of clusters
        -try with M1,M2...Mi and R1, R2,...Ri and see if you can compare the
         variance bwt the maps
    produce
        -SNR as a function number of cluster, one for central bin and one for R500 bin
        -Then do again and change parameters of AR to see what happens 
