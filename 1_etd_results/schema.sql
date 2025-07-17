CREATE TABLE genomes (
                      id INTEGER PRIMARY KEY,
                      genome_name TEXT NOT NULL,
                      organism TEXT);
CREATE TABLE sketch (
                      id INTEGER PRIMARY KEY,
                      sketch_path TEXT NOT NULL);
CREATE TABLE annotations (
                      id INTEGER PRIMARY KEY,
                      genome_id INTEGER REFERENCES genomes(id),
                      gene_name TEXT,
                      amr_annotation TEXT,
                      amr_output_path TEXT,
                      plasmid_annotation TEXT,
                      plasmid_output_path TEXT);
