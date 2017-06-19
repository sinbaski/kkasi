drop procedure if exists prev_row;

delimiter #
create procedure prev_row (
       in today date,
       in tbl varchar(16),
       out day date,
       out high decimal(7,2),
       out low decimal(7, 2),
       out closing decimal(7, 2),
       out volume int
)
begin
	set @m = concat(
	    'select * into @day, @high, @low, @closing, @volume from ',
	    tbl, ' where day < "', today, '" and ',
	    'day >= all (select day from ', tbl,
	    ' where day < "', today, '");'
	    );

	prepare stmt from @m;
	execute stmt;
	deallocate prepare stmt;
	select @day, @high, @low, @closing, @volume into day, high, low, closing, volume;
end#


drop function if exists relret;
delimiter #
create function relret (
       day date,
       tbl varchar(20),
       price decimal(7,2)
) returns decimal
begin
	-- declare day date;
	-- declare tbl varchar(16);
	-- set @day = '2014-09-15';
	-- set @tbl = 'KO_US';
	-- call prev_row(day, tbl, @prev_day,
	-- @high, @low, @closing, @volume);
		
	select (price - @closing)/@closing into @x;
	return @x;
end#

delimiter ;
