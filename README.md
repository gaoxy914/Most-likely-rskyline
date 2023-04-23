# Most-likely-rskyline

## Compliation

Make directory obj and then run makefile.

## Generate synthetic datasets

1. make directories data/ind/ and data/anti/
2. make directories query/weak/ and query/inter/
3. run ./main -genDQ

## Generate nba dataset and car dataset

Run python3 nba.py or python3 car.py then move nba.data or car.dat to data/.

## Run algorithms

Commands format:
```
./main DSA/BBA/LSA data n d alpha query c tau
```
|Parameters|Meaning|Example|
|:---:|---|:---:|
|data|data type|ind, anti, car, nba|
|n|data size|10000, ...|
|d|data dimension|2, 3, 4, ...|
|alpha|probability center|0.2, ...|
|query|query type|weak inter|
|c|number of constraints|1, 2, ...|
|tau|hyper-parameter for LSA|1.2, ...|
