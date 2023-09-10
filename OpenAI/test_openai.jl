println("Testing OpenAI API")

using OpenAI

secret_key = get(ENV,"OPENAI_API_KEY","")
model = "gpt-3.5-turbo"
prompt =  "You are a expert programmer in Fortran, Julia and High Performance computing."

r = create_chat(
    secret_key, 
    model,
    [Dict("role" => "user", "content"=> prompt)]
  )
println(r.response[:choices][begin][:message][:content])
