println("Testing OpenAI API")

using OpenAI

secret_key = get(ENV,"OPENAI_API_KEY","")
model = "gpt-3.5-turbo"
prompt =  "You are going to help me convert Fortran code to Julia."

Msgs = [Dict("role" => "system", "content"=> prompt)]

#Options = Dict("temperature" => "0.1")

r = create_chat(
    secret_key, 
    model,
    Msgs,
    temperature=0.1
  )
println(r.response[:choices][begin][:message][:content])

push!(Msgs, Dict("role" => r.response[:choices][1][:message][:role], "content"=> r.response[:choices][1][:message][:content] ))

prompt = "Can you repeat what I just said?"

push!(Msgs, Dict("role" => "user", "content"=> prompt))

r = create_chat(
    secret_key, 
    model,
    Msgs,
    temperature=0.1
  )
println(r.response[:choices][begin][:message][:content])


# prompt =  "Who are you?"
# r2 = create_chat(
#     secret_key, 
#     model,
#     [Dict("role" => "user", "content"=> prompt)]
#   )
# println(r2.response[:choices][begin][:message][:content])




