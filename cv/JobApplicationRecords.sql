create table jobs (
       title varchar(64),
       company varchar(64),
       recruiter varchar(64),
       salary float,
       city varchar(16),
       country varchar(32),
       link varchar(256),
       cv varchar(64),
       letter varchar(64),
       status enum('viewed', 'applied', 'ToInterview', 'interviewed') default 'viewed',
       primary key (title, company)
);
