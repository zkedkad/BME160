FASTQ = "@EAS139:136:FC706VJ:2:2104:15343:197393"
username = FASTQ.find(":")
password = FASTQ.find(":", (username + 1))
user_id = FASTQ.find(":", (password + 1))
print(FASTQ[1:username])
print(FASTQ[username + 1:password])
print(FASTQ[password + 1:user_id])
