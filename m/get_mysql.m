function mysql = get_mysql()
timefmt = 'yyyy-mm-dd HH:MM:SS';
if exist('/nfs/users3/xiexiaol', 'dir') == 7
    hostname = '83.176.196.41';
    jdbc = ['/nfs/users3/xiexiaol/lib/mysql-connector-java-5.1.26/'...
            'mysql-connector-java-5.1.26-bin.jar'];
    if isempty(strfind(javaclasspath, jdbc))
        javaaddpath(jdbc);
    end
elseif exist('/home/lxb353', 'dir') == 7
    hostname = '83.176.196.41';
    jdbc = ['/home/lxb353/myspace/lib/mysql-connector-java-5.1.33-bin.jar'];
    if isempty(strfind(javaclasspath, jdbc))
        javaaddpath(jdbc);
    end
else
    hostname = 'localhost';
end
mysql = database('avanza', 'sinbaski', 'q1w2e3r4',...
                 'com.mysql.jdbc.Driver', ...
                 sprintf('jdbc:mysql://%s:3306/avanza', hostname));
