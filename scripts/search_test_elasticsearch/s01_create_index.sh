# curl -X PUT "localhost:9200/gnomad?pretty"
curl -X PUT "localhost:9200/gnomad?pretty" -H "Content-Type: application/json" --data-binary @gnomad_mapping.json
curl -X PUT "localhost:9200/gnomad/_settings" -H "Content-Type: application/json" -d '{"index": {"mapping.total_fields.limit": 1000000} }'
curl -X GET "localhost:9200/gnomad?pretty"
