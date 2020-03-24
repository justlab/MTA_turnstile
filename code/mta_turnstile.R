# Get the number of New York City Metropolitan Transit Authority
# (MTA) turnstile entries and exits by station and date.

suppressPackageStartupMessages(
   {library(data.table)
    library(stringr)})

source("../Just_universal/code/download.R")

data.root = "/data-coco/COVID_19"
date.first = as.Date("2020-02-01")

download = function(url, to, f, ...)
    download.update.meta(data.root, url, to, f, ...)

turnstile = function()
  # Field descriptions:
  # https://web.archive.org/web/20190904142924/http://web.mta.info/developers/resources/nyct/turnstile/ts_Field_Description.txt
   {url.root = "http://web.mta.info/developers/data/nyct/turnstile/"
    file.interval.days = 7L

    d = rbindlist(lapply(
        seq(date.first, lubridate::today("America/New_York"),
            by = file.interval.days),
        function(the.date)
           {message("turnstile - ", the.date)
            fname = sprintf("turnstile_%02d%02d%02d.txt",
                year(the.date) %/% 100L, month(the.date), mday(the.date))
            download(
                paste0(url.root, fname),
                file.path("mta_turnstile", fname),
                fread)}))

    message("Processing")

    setnames(d, tolower(colnames(d)))
    setnames(d, "c/a", "ca")

    stopifnot(!anyNA(d))
    stopifnot(all(d[, .(entries, exits)] >= 0))
    d = d[desc != "RECOVR AUD"]
      # Remove recovered-audit rows, which can interrupt otherwise
      # cumulative counts.
    d[, date := lubridate::mdy(date)]
    stopifnot(all(str_detect(d$time, "\\A\\d\\d:\\d\\d:\\d\\d\\z")))
    stopifnot(all(d[,
        by = .(ca, unit, scp, station, date, time), .N]$N == 1))
    setkey(d, ca, unit, scp, station, date, time)

    # The counts of entries and exits are cumulative per (`ca`,
    # `unit`, `scp`, `station`) tuple (at least, I think that's the
    # right tuple), but sometimes the counters are reset. So, we
    # need to subtract each count from the next unless the result
    # would be negative, which indicates a reset. This method can only
    # detect resets imperfectly, but I can't see any other way to do
    # it. Worse, the number to which the counter is reset isn't
    # always 0, so we can't even assume that there have been at least
    # as many entries or exits as the new count.
    #
    # And then there are some cases that are just weird, such as
    # some turnstiles appearing to count backwards for a while,
    # (possibly due to a 64-bit integer being interpreted as 32 bits)
    # and then going back to forwards counting. We'll try to
    # exclude these cases.

    min.obs.n = 10L
    max.neg.diffs.p = .01
    min.pos.diffs.n = 3L
    max.diff = 100e3

    vnames = c("entries", "exits")
    l = lapply(vnames, function(vname)
       {n.dropped = 0
        out = d[, by = .(ca, unit, scp, station),
           {x.old = get(vname)
            diffs = diff(x.old)
            if (.N >= min.obs.n &&
                    mean(diffs < 0) <= max.neg.diffs.p &&
                    sum(diffs > 0) >= min.pos.diffs.n &&
                    all(diffs <= max.diff))
                .(date = date[-1], count = diffs)
            else
               {n.dropped <<- n.dropped + 1
                NULL}}]
        message(sprintf("%s: dropped %d tuples", vname, n.dropped))
        out = out[count >= 0]

        out[, by = .(station, date), .(
            units = length(unique(paste(ca, unit, scp))),
            count = sum(count))]})

    merge(l[[1]], l[[2]], all = T,
        by = c("station", "date"),
        suffixes = paste0(".", vnames))}
