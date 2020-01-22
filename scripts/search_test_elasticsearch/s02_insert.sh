# curl -X POST "localhost:9200/gnomad/_bulk?refresh" -H "Content-Type: application/json" --data-binary @$1 >> aa
curl -X POST "127.0.0.1:9200/gnomad/_bulk?refresh" -H "Content-Type: application/json" --data-binary @$1 >> aa
# curl -X POST "compute-e-16-233:9200/gnomad/_bulk?refresh" -H "Content-Type: application/json" --data-binary @$1 >> aa
# curl -X POST "10.120.16.233:9200/gnomad/_bulk?refresh" -H "Content-Type: application/json" --data-binary @$1 >> aa
# curl -X POST "login06.o2.rc.hms.harvard.edu:9200/gnomad/_bulk?refresh" -H "Content-Type: application/json" --data-binary @$1 
# curl -X POST "134.174.159.26:9200/gnomad/_bulk?refresh" -H "Content-Type: application/json" --data-binary @$1 
# curl -X POST "compute-e-16-233.o2.rc.hms.harvard.edu:9200/gnomad/_bulk?refresh" -H "Content-Type: application/json" --data-binary @$1 >> aa
