SELECT doc_id, count(distinct molregno) cnt into temporary table doc_compound_counts \
      from docs join compound_records using (doc_id) join compound_structures using (molregno) \
      group by (doc_id);

SELECT count(*) from doc_compound_counts where cnt>=20 and cnt<=100;
