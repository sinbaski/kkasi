function mysql = get_mysql()
timefmt = 'yyyy-mm-dd HH:MM:SS';
if exist('/nfs/users3/xiexiaol', 'dir') == 7
    hostname = '85.228.154.211';
    jdbc = ['/nfs/users3/xiexiaol/lib/mysql-connector-java-5.1.26/'...
            'mysql-connector-java-5.1.26-bin.jar'];
elseif exist('/home/lxb353', 'dir') == 7
    hostname = '85.228.154.211';
    jdbc = ['/home/lxb353/myspace/lib/mysql-connector-java-5.1.33-bin.jar'];
else
    hostname = 'localhost';
    jdbc = ['/home/xxie/.matlab/mysql-connector-java-5.1.35-bin.jar'];
end

if isempty(strfind(javaclasspath, jdbc))
    javaaddpath(jdbc);
end
mysql = database('avanza', 'sinbaski', 'q1w2e3r4',...
                 'com.mysql.jdbc.Driver', ...
                 sprintf('jdbc:mysql://%s:3306/avanza', hostname));
