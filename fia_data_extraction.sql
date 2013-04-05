-- Species abundance distribution data from Forest Inventory
-- and Analysis data.
-- Builds on https://github.com/weecology/white-etal-2012-ecology/blob/master/mete_sads_data.py
-- adding additional filtering based on human disturbance and
-- evidence of cutting

DROP TABLE IF EXISTS queries.fia_survey1;
DROP TABLE IF EXISTS queries.fia_plot1;
DROP TABLE IF EXISTS queries.fia_plot2;
DROP TABLE IF EXISTS queries.fia_cond1;
DROP TABLE IF EXISTS queries.fia_plot3;
DROP TABLE IF EXISTS queries.fia_plot4;
DROP TABLE IF EXISTS queries.fia_tree1;
DROP TABLE IF EXISTS queries.fia_tree2;
DROP TABLE IF EXISTS queries.fia_tree3;

CREATE TABLE queries.fia_survey1
SELECT SURVEY.CN
FROM FIA.SURVEY
WHERE SURVEY.ANN_INVENTORY = '"Y"';

CREATE TABLE queries.fia_plot1
SELECT PLOT.CN, PLOT.SRV_CN, PLOT.CTY_CN, PLOT.INVYR,
PLOT.STATECD, PLOT.UNITCD, PLOT.COUNTYCD, PLOT.PLOT, PLOT.LAT,
PLOT.LON, PLOT.ELEV
FROM FIA.PLOT
WHERE PLOT.PLOT_STATUS_CD = 1 AND PLOT.KINDCD > 0 AND
PLOT.KINDCD < 4 AND
PLOT.DESIGNCD IN (1, 311, 312, 313, 314, 328, 220, 240) AND
PLOT.MANUAL >= 1 AND PLOT.QA_STATUS = 1 AND
PLOT.SAMP_METHOD_CD = 1 AND PLOT.INVYR < 3000;

CREATE TABLE queries.fia_plot2
SELECT fia_plot1.CN, fia_plot1.CTY_CN, fia_plot1.INVYR,
fia_plot1.STATECD, fia_plot1.UNITCD, fia_plot1.COUNTYCD,
fia_plot1.PLOT, fia_plot1.LAT, fia_plot1.LON, fia_plot1.ELEV
FROM queries.fia_plot1 INNER JOIN queries.fia_survey1
ON fia_plot1.SRV_CN = fia_survey1.CN
ORDER BY fia_plot1.STATECD, fia_plot1.UNITCD,
fia_plot1.COUNTYCD, fia_plot1.PLOT, fia_plot1.INVYR DESC;

CREATE TABLE queries.fia_cond1
SELECT COND.PLT_CN
FROM FIA.COND
WHERE ((COND.STDORGCD)="0" OR (COND.STDORGCD)='')
AND ((COND.TRTCD1) = '' OR (COND.TRTCD1)="0")
AND ((COND.TRTCD2) = '' OR (COND.TRTCD2)="0")
AND ((COND.TRTCD3) = '' OR (COND.TRTCD3)="0")
AND (COND.DSTRBCD1 <> 80)
AND (COND.DSTRBCD2 <> 80)
AND (COND.DSTRBCD3 <> 80)
AND ((COND.STUMP_CD_PNWRS = '"N"') OR (COND.STUMP_CD_PNWRS = '""')) -- the quotes themselves got imported into the database; sigh
GROUP BY COND.PLT_CN;

CREATE TABLE queries.fia_plot3
SELECT fia_plot2.CN, fia_plot2.STATECD, fia_plot2.UNITCD,
fia_plot2.COUNTYCD, fia_plot2.PLOT,
MAX(fia_plot2.INVYR) AS INVYR, AVG(fia_plot2.LAT) AS LAT,
AVG(fia_plot2.LON) AS LON, AVG(fia_plot2.ELEV) AS ELEV
FROM queries.fia_plot2
GROUP BY fia_plot2.STATECD, fia_plot2.UNITCD,
fia_plot2.COUNTYCD, fia_plot2.PLOT;

CREATE TABLE queries.fia_plot4
SELECT fia_plot3.*
FROM queries.fia_plot3
INNER JOIN queries.fia_cond1 ON fia_plot3.CN = fia_cond1.PLT_CN;

CREATE TABLE queries.fia_tree1
SELECT TREE.PLT_CN, TREE.STATECD, TREE.UNITCD, TREE.COUNTYCD,
TREE.PLOT, TREE.SPCD
FROM FIA.TREE
WHERE TREE.STATUSCD = 1;

CREATE TABLE queries.fia_tree2
SELECT ((fia_tree1.STATECD*10000000000) +
(fia_tree1.UNITCD*1000000000) + (fia_tree1.COUNTYCD*1000000) +
fia_tree1.PLOT) AS PlotID, fia_tree1.PLT_CN, fia_tree1.SPCD
FROM queries.fia_tree1 INNER JOIN queries.fia_plot4
ON fia_tree1.PLT_CN = fia_plot4.CN
WHERE fia_plot4.INVYR < 3000;

CREATE TABLE queries.fia_tree3
SELECT fia_tree2.PlotID, fia_tree2.PLT_CN, fia_tree2.SPCD,
COUNT(fia_tree2.SPCD) AS AB
FROM queries.fia_tree2
GROUP BY fia_tree2.PlotID, fia_tree2.PLT_CN, fia_tree2.SPCD;

SELECT * FROM queries.fia_tree3
INTO OUTFILE '/tmp/fia_spab.csv'
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n';
