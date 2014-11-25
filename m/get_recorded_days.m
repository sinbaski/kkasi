function days = get_recorded_days(company, first_day, last_day)
mysql = get_mysql();
stmt = sprintf(['select distinct(date(tid)) from %s where date(tid) between' ...
                ' "%s" and "%s" order by tid;'], company, ...
               first_day, last_day);
days = fetch(mysql, stmt);

close(mysql);
