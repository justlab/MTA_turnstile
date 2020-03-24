suppressPackageStartupMessages(
   {library(data.table)
    library(stringr)
    library(pbapply)
    library(sf)})

source("../Just_universal/code/pairmemo.R")
source("../Just_universal/code/download.R")

data.root = "/data-coco/COVID_19"
pairmemo.dir = file.path(data.root, "pairmemo")
date.first = as.Date("2014-12-27")

crs.lonlat = 4326 # https://epsg.io/4326
crs.us.atlas = 2163 # https://epsg.io/2163

comma = scales::comma

download = function(url, to, f, ...)
    download.update.meta(data.root, url, to, f, ...)

varn = function(v, na.rm = F)
# Variance with n instead of (n - 1) in the denominator.
    mean((v - mean(v, na.rm = na.rm))^2, na.rm = na.rm)
sdn = function(v, na.rm = F)
# Standard deviation with n instead of (n - 1) in the denominator.
    sqrt(varn(v, na.rm))

turnstile = function()
  # Get the number of New York City Metropolitan Transit Authority
  # (MTA) turnstile entries and exits by station and date.
  # Field descriptions:
  # https://web.archive.org/web/20190904142924/http://web.mta.info/developers/resources/nyct/turnstile/ts_Field_Description.txt
   {url.root = "http://web.mta.info/developers/data/nyct/turnstile/"
    file.interval.days = 7L

    message("Reading turnstile files")
    d = rbindlist(pblapply(
        seq(date.first, lubridate::today("America/New_York"),
            by = file.interval.days),
        function(the.date)
           {fname = sprintf("turnstile_%02d%02d%02d.txt",
              # `strftime` is platform-dependent, so it's a bit
              # safer to do this by hand.
                year(the.date) %% 1000L,
                month(the.date),
                mday(the.date))
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
    # `unit`, `scp`) tuple (at least, I think that's the
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

    # First, remove some special cases, where counts jump back
    # and forth between observations. Probably two turnstiles
    # were being counted as one.

    d = d[!(
        (ca == "PTH12" & unit == "R542" & scp == "00-04-00" &
           date %in% (as.Date("2019-04-27") + (0:3))) |
        (ca == "PTH13" & unit == "R541" & scp == "00-00-04" &
           date %in% (as.Date("2016-02-02") + (0:2))))]

    # Now handle the rest.

    min.obs.n = 10L
    max.neg.diffs.p = .01
    min.pos.diffs.n = 3L
    max.diff = 100e3

    vnames = c("entries", "exits")
    l = lapply(vnames, function(vname)
       {message("Decumulating - ", vname)
        n.sources = nrow(d[, by = .(ca, unit, scp, station), 1])
        n.sources.dropped = 0
        n.obs = nrow(d)
        n.obs.dropped = 0
        bar = txtProgressBar(style = 3, min = 0, max = n.sources)
        out = d[, by = .(ca, unit, scp),
           {setTxtProgressBar(bar, .GRP)
            x.old = get(vname)
            diffs = diff(x.old)
            if (.N >= min.obs.n &&
                    mean(diffs < 0) <= max.neg.diffs.p &&
                    sum(diffs > 0) >= min.pos.diffs.n &&
                    all(diffs <= max.diff))
               {date.diffs = as.integer(diff(date))
                i.keep = date.diffs <= 1 & diffs >= 0
                n.obs.dropped <<- n.obs.dropped + sum(!i.keep)
                .(date = date[-1][i.keep], count = diffs[i.keep])}
            else
               {n.sources.dropped <<- n.sources.dropped + 1
                n.obs.dropped <<- n.obs.dropped + .N
                NULL}}]
        close(bar)
        message(sprintf("Dropped %s sources of %s",
            comma(n.sources.dropped), comma(n.sources)))
        message(sprintf("Dropped %s observations of %s",
            comma(n.obs.dropped), comma(n.obs)))

        message("Summarizing")
        out[, by = .(ca, date), .(
            sources = length(unique(paste(unit, scp))),
            count = sum(count))]})

    message("Merging")
    counts = merge(l[[1]], l[[2]], all = T,
        by = c("ca", "date"),
        suffixes = paste0(".", vnames))

    # Observations with the same `ca` but different station names
    # should be referring to the same station. Use the latest name for
    # each station. Beware that different `ca`s can have the same
    # station name.
    message("Setting station names")
    stations = d[, by = ca,
        .(station.name = station[which.max(date)])]

    # Geolocate stations.
    locations = download(
        "https://raw.githubusercontent.com/chriswhong/nycturnstiles/master/geocoded.csv",
        "mta_turnstile_locations.csv",
        fread, header = F, col.names = c(
            "unit", "ca", "station", "linename", "division",
            "lat", "lon"))
    setkey(locations, ca)
    stations[, c("lon", "lat") :=
        locations[.(stations$ca), .(lon, lat)]]
    # Drop unlocated stations and their associated observations.
    message(sprintf("Dropping %d unlocated stations ",
        length(is.na(stations$lon))))
    stations = stations[!is.na(lon)]
    message(sprintf("Dropping %s observations",
        counts[, comma(sum(!(ca %in% stations$ca)))]))
    counts = counts[ca %in% stations$ca]

    list(counts = counts, stations = stations)}
turnstile = pairmemo(turnstile, pairmemo.dir)

station.boros = function()
  # Return a factor of boro names, with one boro for each station.
   {boros = download(
        "https://data.cityofnewyork.us/api/geospatial/tqmj-j8zm?method=export&format=Original",
        "nyc_boros_shapefile.zip",
        function(p) read_sf(paste0("/vsizip/", p, "/nybb_20a")))
    stations = st_as_sf(turnstile()$stations,
        coords = c("lon", "lat"), crs = crs.lonlat)
    result = do.call(st_intersects, lapply(list(stations, boros),
        function(x) st_transform(x, crs = crs.us.atlas)))
    factor(sapply(result, function(x)
        if (length(x) == 0)
            NA
        else if (length(x) == 1)
            boros[x,]$BoroName
        else
            stop()))}
station.boros = pairmemo(station.boros, pairmemo.dir)

relative.subway.usage = function(the.year)
  # For each date in the given year, compute a z-score of subway usage
  # comparing that date to all other days that occur on the same month
  # and day of the week, but in different years. The results are given
  # per boro and also pooled across boros.
   {counts = turnstile()$counts
    counts[, boro := station.boros()[match(ca, turnstile()$stations$ca)]]
    # Drop observations that we can't assign a boro to.
    counts = counts[!is.na(boro)]

    d = rbindlist(lapply(c(F, T), function(boro.split)
       {d = copy(counts)
        if (!boro.split)
            d[, boro := "all"]
        d = d[, by = .(date, boro), .(uses =
            sum(count.entries, na.rm = T) +
            sum(count.exits, na.rm = T))]
        typical = d[
            year(date) != the.year &
                date < lubridate::make_date(year(date), 4, 1),
            keyby = .(boro, m = month(date), w = wday(date)),
            .(mean.t = mean(uses), sd.t = sdn(uses))]
        d[year(date) == the.year, .(date, boro, usage =
           {x = data.table(boro, month(date), wday(date))
            x = typical[.(x)]
            (uses - x$mean.t) / x$sd.t})]}))

    d = dcast(d, date ~ boro, value.var = "usage")
    setnames(d, make.names(colnames(d)))

    d}
