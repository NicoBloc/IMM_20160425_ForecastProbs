library(decon)

dRaw <- read.xls("../../Data/toyData.xlsx", sheet=1)[, 1]

sig=1.5

y2 <- dRaw[!(dRaw %in% c(6, 40))]
dd <- DeconPdf(y=y2, sig=sig, error='normal', from=min(y2), to=max(y2), n=(diff(range(y2)) + 1), bw=0.75)

h <- hist(y2, breaks=(min(y2) - 0.5):(max(y2) + 0.5), freq=FALSE)
plot(h$mids, h$density, type='o')
points(dd$x, dd$y, col=4, type='o')

bw <- 5
err <- dnorm(-bw:bw, sd=sig)
err <- err / sum(err)
err
nErr <- 2 * bw + 1
yConv <- convolve(dd$y, err, type='open')


points((min(dd$x) - bw):(max(dd$x) + bw), yConv, col=2, type='o')
legend('topleft', legend=c('observed diameter', 'estimated true diameter', 'reconvoluted'), col=c(1, 4, 2), lty=1)




# With deamer

dd <- deamerKE(y=y2, mu=0, sigma=sig, noise.type='Gaussian', from=min(y2), to=max(y2), grid.length=(diff(range(y2)) + 1))
s